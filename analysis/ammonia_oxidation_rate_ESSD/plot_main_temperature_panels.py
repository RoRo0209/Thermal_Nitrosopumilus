#!/usr/bin/env python3
"""Create compact, panel-ready temperature-response figures for the main text."""

# Release path adapter: data and outputs live in a separate user workspace.
import os as _release_os
from pathlib import Path as _ReleasePath
if not _release_os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the analysis workspace")
_release_root = _ReleasePath(_release_os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
_release_os.chdir(_release_root / 'ammonia_oxidation_rate_ESSD')


import os
from pathlib import Path

HERE = _release_root / "ammonia_oxidation_rate_ESSD"
DATA = HERE / "Database of nitrification and nitrifiers in the global ocean.xlsx"
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".matplotlib"))

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

mpl.use("Agg")
import matplotlib.pyplot as plt

AOA_TEAL = "#167C80"
AOA_POINT = "#9AAAB5"
BLUE = "#2F7FAE"
ORANGE = "#D9653B"
BLACK = "#1A1A1A"
GREY = "#606A70"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def style_axis(ax):
    for spine in ax.spines.values():
        spine.set_color(BLACK)
        spine.set_linewidth(0.7)
    ax.tick_params(direction="out", colors=BLACK, pad=2)
    ax.set_facecolor("white")


def save_panel(fig, stem):
    for suffix, kwargs in (
        ("pdf", {}),
        ("svg", {}),
        ("png", {"dpi": 600}),
        ("tiff", {"dpi": 600}),
    ):
        fig.savefig(HERE / f"{stem}.{suffix}", bbox_inches="tight", **kwargs)


def binned_summary(frame, value_col, min_n):
    bins = np.arange(-2, 32, 2)
    rows = []
    for left, right in zip(bins[:-1], bins[1:]):
        values = frame.loc[
            (frame["temp"] >= left) & (frame["temp"] < right), value_col
        ]
        if len(values) >= min_n:
            rows.append({
                "temp": (left + right) / 2,
                "mean": values.mean(),
                "se": values.sem(),
                "n": len(values),
            })
    return pd.DataFrame(rows)


def load_data():
    aoa = pd.read_excel(DATA, sheet_name="ammonia oxidizer")
    aoa.columns = [str(c).strip() for c in aoa.columns]
    aoa = aoa.rename(columns={
        "Depth (m)": "depth",
        "total AOA (copies L-1)": "total_aoa",
        "shallow clade AOA (copies L-1)": "shallow_aoa",
        "Temperature (°C)": "temp",
    })
    aoa["aoa"] = aoa["total_aoa"].fillna(aoa["shallow_aoa"])
    for col in ("depth", "temp", "aoa"):
        aoa[col] = pd.to_numeric(aoa[col], errors="coerce")
    aoa = aoa[
        (aoa["depth"] <= 200)
        & aoa["temp"].notna()
        & aoa["aoa"].notna()
        & (aoa["aoa"] > 0)
    ].copy()
    aoa["log10_aoa"] = np.log10(aoa["aoa"])

    rate = pd.read_excel(DATA, sheet_name="volumetric ammonia oxidation")
    rate.columns = [str(c).strip() for c in rate.columns]
    rate = rate.rename(columns={
        "Depth (m)": "depth",
        "Rate (nmol N L-1 d-1)": "rate",
        "Temperature (°C)": "temp",
        "Ammonium+ammonia (μM)": "nh4",
    })
    for col in ("depth", "temp", "rate", "nh4"):
        rate[col] = pd.to_numeric(rate[col], errors="coerce")
    rate = rate[
        (rate["depth"] <= 200)
        & rate["temp"].notna()
        & rate["rate"].notna()
        & (rate["rate"] > 0)
        & rate["nh4"].notna()
    ].copy()
    rate["log10_rate"] = np.log10(rate["rate"])
    threshold = rate["nh4"].median()
    rate["nh4_group"] = np.where(
        rate["nh4"] <= threshold, "Low NH4", "High NH4"
    )
    return aoa, rate, threshold


def plot_aoa_panel(aoa):
    summary = binned_summary(aoa, "log10_aoa", min_n=5)
    fig, ax = plt.subplots(figsize=(65 / 25.4, 56 / 25.4))

    ax.scatter(
        aoa["temp"], aoa["log10_aoa"],
        s=5.0, color=AOA_POINT, alpha=0.24, edgecolors="none", rasterized=True,
    )
    ax.errorbar(
        summary["temp"], summary["mean"], yerr=summary["se"],
        color=AOA_TEAL, linewidth=1.55, marker="o", markersize=3.4,
        markerfacecolor="white", markeredgecolor=AOA_TEAL, markeredgewidth=0.85,
        capsize=1.8, elinewidth=0.8, zorder=4,
        label="Binned mean ± s.e.",
    )

    ax.set_xlim(-3, 33)
    ax.set_ylim(2, 9.1)
    ax.set_xticks([0, 10, 20, 30])
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel(r"log$_{10}$(AOA amoA copies L$^{-1}$)")
    ax.legend(
        loc="upper right", fontsize=6.0, handlelength=1.7,
        borderaxespad=0.35, labelcolor=AOA_TEAL,
    )
    ax.text(
        0.04, 0.05, f"n = {len(aoa):,}", transform=ax.transAxes,
        fontsize=6.2, color=GREY,
    )
    style_axis(ax)
    fig.subplots_adjust(left=0.22, right=0.98, bottom=0.20, top=0.97)
    save_panel(fig, "Fig_main_panel_AOA_abundance_vs_temp")
    plt.close(fig)

    aoa[["temp", "aoa", "log10_aoa"]].to_csv(
        HERE / "SourceData_main_panel_AOA_abundance_vs_temp.csv", index=False
    )


def format_p(p):
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}"


def plot_rate_panel(rate, threshold):
    low = rate[rate["nh4"] <= threshold].copy()
    high = rate[rate["nh4"] > threshold].copy()
    low_summary = binned_summary(low, "log10_rate", min_n=3)
    high_summary = binned_summary(high, "log10_rate", min_n=3)
    rho_low, p_low = stats.spearmanr(low["temp"], low["log10_rate"])
    rho_high, p_high = stats.spearmanr(high["temp"], high["log10_rate"])

    fig, ax = plt.subplots(figsize=(65 / 25.4, 56 / 25.4))
    ax.scatter(
        high["temp"], high["log10_rate"], s=6, color=ORANGE,
        alpha=0.18, edgecolors="none", rasterized=True,
    )
    ax.scatter(
        low["temp"], low["log10_rate"], s=6, color=BLUE,
        alpha=0.20, edgecolors="none", rasterized=True,
    )

    for summary, color in ((high_summary, ORANGE), (low_summary, BLUE)):
        ax.errorbar(
            summary["temp"], summary["mean"], yerr=summary["se"],
            color=color, linewidth=1.55, marker="o", markersize=3.4,
            markerfacecolor="white", markeredgecolor=color, markeredgewidth=0.85,
            capsize=1.8, elinewidth=0.8, zorder=4,
        )

    ax.set_xlim(-3, 33)
    ax.set_ylim(-4.3, 3.6)
    ax.set_xticks([0, 10, 20, 30])
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel(r"log$_{10}$(AO rate, nmol N L$^{-1}$ d$^{-1}$)")

    handles = [
        Line2D([0], [0], color=BLUE, marker="o", markerfacecolor="white",
               markeredgecolor=BLUE, linewidth=1.5, markersize=3.2),
        Line2D([0], [0], color=ORANGE, marker="o", markerfacecolor="white",
               markeredgecolor=ORANGE, linewidth=1.5, markersize=3.2),
    ]
    labels = [
        rf"Low NH$_4^+$ ≤{threshold:.2f} μM",
        rf"High NH$_4^+$ >{threshold:.2f} μM",
    ]
    ax.legend(
        handles, labels, loc="upper center", ncol=2, fontsize=5.7,
        handlelength=1.4, columnspacing=0.8, handletextpad=0.35,
        borderaxespad=0.35,
    )
    style_axis(ax)
    fig.subplots_adjust(left=0.22, right=0.98, bottom=0.20, top=0.97)
    save_panel(fig, "Fig_main_panel_AO_rate_vs_temp_by_NH4")
    plt.close(fig)

    rate[["temp", "rate", "log10_rate", "nh4", "nh4_group"]].to_csv(
        HERE / "SourceData_main_panel_AO_rate_vs_temp_by_NH4.csv", index=False
    )
    return {
        "nh4_threshold": threshold,
        "low_n": len(low),
        "low_rho": rho_low,
        "low_p": p_low,
        "high_n": len(high),
        "high_rho": rho_high,
        "high_p": p_high,
    }


def make_preview():
    """Rasterize both Python-generated PDFs into a side-by-side QA preview."""
    import fitz

    names = [
        "Fig_main_panel_AOA_abundance_vs_temp.pdf",
        "Fig_main_panel_AO_rate_vs_temp_by_NH4.pdf",
    ]
    images = []
    for name in names:
        doc = fitz.open(HERE / name)
        pix = doc[0].get_pixmap(matrix=fitz.Matrix(2.4, 2.4), alpha=False)
        arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)
        images.append(arr)

    fig, axes = plt.subplots(1, 2, figsize=(150 / 25.4, 62 / 25.4))
    for ax, image, letter in zip(axes, images, ("a", "b")):
        ax.imshow(image)
        ax.axis("off")
        ax.text(0.005, 0.995, letter, transform=ax.transAxes, ha="left", va="top",
                fontsize=8, fontweight="bold", color=BLACK)
    fig.subplots_adjust(left=0.005, right=0.995, bottom=0.01, top=0.99, wspace=0.08)
    fig.savefig(HERE / "Fig_main_temperature_panels_preview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    aoa, rate, threshold = load_data()
    plot_aoa_panel(aoa)
    stats_out = plot_rate_panel(rate, threshold)
    make_preview()
    print(f"AOA n = {len(aoa)}")
    for key, value in stats_out.items():
        print(f"{key} = {value}")


if __name__ == "__main__":
    main()
