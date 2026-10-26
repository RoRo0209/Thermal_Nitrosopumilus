from pathlib import Path
from collections import Counter
import re
import zlib

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch, Rectangle
from PIL import Image


import os
if not os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to a separate workspace containing the documented input tables")
ROOT = Path(os.environ["THERMAL_PROJECT_ROOT"]).resolve()
OUT = ROOT / "ALE_2" / "Fig5_asymmetric_accessory"
FAMILY_MATRIX = ROOT / "ALE_2" / "genome_x_family_copy_number.tsv"
KO_MATRIX = ROOT / "ALE_2" / "full_ko_analysis" / "genome_by_KO_copy_number.tsv"
MAPPING = ROOT / "ALE_2" / "tree-gen-name-gff.xlsx"
METADATA = ROOT / "genome_character" / "functional_gene_counts_all_np_processed.csv"

GROUPS = ["highT", "lowT", "other_NP"]
GROUP_LABELS = {"highT": "HighT", "lowT": "LowT", "other_NP": "Other"}
GROUP_COLORS = {"highT": "#D97732", "lowT": "#4F93B5", "other_NP": "#9A9A9A"}
GROUP_BOX_COLORS = {
    "highT": "#F1CBB1",
    "lowT": "#BCD6E2",
    "other_NP": "#D9D9D9",
}
BACKGROUND_COLOR = "#FFFFFF"
TEXT_COLOR = "#000000"
PREVALENCE_COLORS = ["#F7F5F8", "#D8CBE1", "#8A6AA0", "#3F2454"]
AUDIT_COLORS = {
    **GROUP_COLORS,
    **{f"{key}_box": value for key, value in GROUP_BOX_COLORS.items()},
    "prevalence_accent": "#3F2454",
}
PREVALENCE_CMAP = LinearSegmentedColormap.from_list(
    "prevalence", PREVALENCE_COLORS
)


def prevalence_hex(value):
    return mpl.colors.to_hex(
        PREVALENCE_CMAP(float(value) / 100), keep_alpha=False
    ).upper()

PANEL_A = [
    ("glnBI", r"$\it{glnB}$"),
    ("cysO", r"$\it{cysK}$"),
    ("betB", r"$\it{gabD}$"),
    ("group_1296", r"$\it{fabG}$"),
    ("cobB1", r"$\it{cobB}$"),
    ("cpdA", r"$\it{cpdA}$"),
    ("group_1327", "UDG family"),
]
PANEL_B = [("group_1948", r"$\it{spoT}$/" + "\nMESH1-like")]
PANEL_C = [("K13993", r"$\it{HSP20}$-like"), ("K00952", r"$\it{nadM}$")]
PANEL_A_GROUPS = [
    (0, 1, "Nitrogen / amino-acid\nmetabolism"),
    (2, 3, "Carbon & energy\nmetabolism"),
    (4, 5, "Regulation /\nsignaling"),
    (6, 6, "Genome\nmaintenance"),
]


def load_matrix(path):
    matrix = pd.read_csv(path, sep="\t")
    mapping = (
        pd.read_excel(MAPPING, sheet_name="Sheet1")[["gff_name", "GCA"]]
        .dropna()
        .drop_duplicates()
    )
    metadata = (
        pd.read_csv(METADATA)[["Genome_ID", "Clade"]]
        .drop_duplicates()
    )
    merged = (
        matrix.merge(mapping, left_on="genome", right_on="gff_name", how="left")
        .merge(metadata, left_on="GCA", right_on="Genome_ID", how="inner")
    )
    observed = merged["Clade"].value_counts().to_dict()
    expected = {"highT": 11, "lowT": 23, "other_NP": 74}
    if observed != expected or merged["Genome_ID"].nunique() != 108:
        raise RuntimeError(f"Genome-to-clade mapping failed: {observed}")
    return merged


def prevalence_data(matrix, definitions, panel):
    rows = []
    for family, label in definitions:
        if family not in matrix.columns:
            raise KeyError(f"Missing family in raw matrix: {family}")
        for group in GROUPS:
            values = matrix.loc[matrix["Clade"] == group, family].astype(float)
            present = int((values > 0).sum())
            rows.append({
                "panel": panel,
                "feature_label": label.replace(r"$\it{", "").replace("}$", "")
                    .replace("\n", " "),
                "feature_id": family,
                "group": GROUP_LABELS[group],
                "genome_id": "",
                "measure": "prevalence_percent",
                "value": 100 * present / len(values),
                "present_genomes": present,
                "total_genomes": len(values),
            })
    return pd.DataFrame(rows)


def copy_number_data(matrix):
    rows = []
    for ko, label in PANEL_C:
        if ko not in matrix.columns:
            raise KeyError(f"Missing KO in copy-number matrix: {ko}")
        clean_label = label.replace(r"$\it{", "").replace("}$", "")
        for _, row in matrix.iterrows():
            rows.append({
                "panel": "c",
                "feature_label": clean_label,
                "feature_id": ko,
                "group": GROUP_LABELS[row["Clade"]],
                "genome_id": row["Genome_ID"],
                "measure": "copy_number_per_genome",
                "value": float(row[ko]),
                "present_genomes": np.nan,
                "total_genomes": np.nan,
            })
    return pd.DataFrame(rows)


def draw_prevalence(
    ax, data, definitions, title, panel_letter, highlight=None,
    functional_groups=None,
):
    labels = [label for _, label in definitions]
    groups = [GROUP_LABELS[g] for g in GROUPS]
    values = np.array([
        [
            data.loc[
                (data["feature_id"] == family) & (data["group"] == group),
                "value",
            ].iloc[0]
            for family, _ in definitions
        ]
        for group in groups
    ])
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            ax.add_patch(Rectangle(
                (col - 0.5, row - 0.5), 1, 1,
                facecolor=prevalence_hex(values[row, col]),
                edgecolor="#FFFFFF", linewidth=0.25,
            ))
            color = "#FFFFFF" if values[row, col] >= 58 else "#222222"
            ax.text(
                col, row, f"{values[row, col]:.0f}", ha="center", va="center",
                fontsize=6.6, color=color,
            )
    ax.set_xlim(-0.5, len(labels) - 0.5)
    ax.set_ylim(len(groups) - 0.5, -0.5)
    ax.set_xticks(range(len(labels)), labels, fontsize=7.3)
    if len(labels) > 3:
        plt.setp(ax.get_xticklabels(), rotation=32, ha="right",
                 rotation_mode="anchor")
    ax.set_yticks(range(len(groups)), groups, fontsize=7.4)
    ax.tick_params(axis="both", length=0, pad=3)
    ax.text(
        0, 1.23, title, transform=ax.transAxes, fontsize=8.2,
        fontweight="bold", ha="left", va="top",
    )
    letter_x = -0.30 if len(labels) == 1 else -0.10
    ax.text(
        letter_x, 1.23, panel_letter, transform=ax.transAxes,
        fontsize=10, fontweight="bold", ha="left", va="top",
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    if highlight is not None:
        row, col = highlight
        ax.add_patch(Rectangle(
            (col - 0.49, row - 0.49), 0.98, 0.98,
            fill=False, lw=1.5, edgecolor="#3F2454",
        ))
    if functional_groups:
        n_features = len(labels)
        for start, end, group_label in functional_groups:
            left = (start + 0.08) / n_features
            right = (end + 0.92) / n_features
            midpoint = (left + right) / 2
            ax.plot(
                [left, right], [1.035, 1.035], transform=ax.transAxes,
                color="#A7A7A7", lw=1.0, solid_capstyle="butt",
                clip_on=False,
            )
            ax.text(
                midpoint, 1.055, group_label, transform=ax.transAxes,
                ha="center", va="bottom", fontsize=5.2, color="#555555",
                linespacing=0.9, clip_on=False,
            )
    return mpl.cm.ScalarMappable(
        norm=Normalize(0, 100), cmap=PREVALENCE_CMAP
    )


def draw_copy_number(ax, data):
    rng = np.random.default_rng(17)
    group_offsets = {"HighT": -0.24, "LowT": 0.0, "Other": 0.24}
    for feature_index, (ko, feature_label) in enumerate(PANEL_C):
        clean_label = feature_label.replace(r"$\it{", "").replace("}$", "")
        for group in ["HighT", "LowT", "Other"]:
            values = data.loc[
                (data["feature_id"] == ko) & (data["group"] == group), "value"
            ].to_numpy()
            position = feature_index + group_offsets[group]
            box = ax.boxplot(
                values, positions=[position], widths=0.17, patch_artist=True,
                showfliers=False, whis=1.5,
                medianprops={"color": "#1F1F1F", "linewidth": 1.0},
                whiskerprops={"color": "#555555", "linewidth": 0.8},
                capprops={"color": "#555555", "linewidth": 0.8},
                boxprops={"edgecolor": "#555555", "linewidth": 0.8},
            )
            box["boxes"][0].set_facecolor(GROUP_BOX_COLORS[
                {"HighT": "highT", "LowT": "lowT", "Other": "other_NP"}[group]
            ])
            box["boxes"][0].set_alpha(1.0)
            jitter = rng.uniform(-0.055, 0.055, len(values))
            ax.scatter(
                np.full(len(values), position) + jitter, values, s=9,
                facecolor=GROUP_COLORS[
                    {"HighT": "highT", "LowT": "lowT", "Other": "other_NP"}[group]
                ],
                edgecolor="#FFFFFF", linewidth=0.25, alpha=1.0, zorder=3,
            )
    ax.set_xticks(range(len(PANEL_C)), [label for _, label in PANEL_C], fontsize=7.6)
    ax.set_ylabel("Copy number per genome", fontsize=7.6)
    ax.text(
        0, 1.23, "Copy-number differences in\nretained functions",
        transform=ax.transAxes, fontsize=8.2, fontweight="bold",
        ha="left", va="top",
    )
    ax.text(
        -0.13, 1.23, "c", transform=ax.transAxes, fontsize=10,
        fontweight="bold", ha="left", va="top",
    )
    ax.set_xlim(-0.52, len(PANEL_C) - 0.48)
    ymax = int(np.ceil(data["value"].max()))
    ax.set_ylim(-0.12, max(2.4, ymax + 0.4))
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.65)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(0.8)
    ax.tick_params(labelsize=7, width=0.7, length=3)
    handles = [
        Patch(facecolor=GROUP_COLORS[g], edgecolor="#FFFFFF", linewidth=0,
              alpha=1.0,
              label=GROUP_LABELS[g])
        for g in GROUPS
    ]
    ax.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.23),
        ncol=3, frameon=False, fontsize=6.8, handlelength=1.1,
        handletextpad=0.35, columnspacing=0.8,
    )


def write_legend():
    text = """# Figure 5 legend draft

**Fig. 5 | Asymmetric accessory-gene repertoires associated with thermal differentiation.**
**a,** Prevalence of representative gene families that are retained in HighT and
other *Nitrosopumilus* genomes but depleted in LowT, grouped by broad functional
roles. Values are percentages of genomes containing at least one family member.
**b,** The structure-supported
*spoT*/MESH1-like family is concentrated in LowT genomes and is nearly absent
from HighT and other *Nitrosopumilus*. **c,** Total copy number per genome for
*HSP20*-like (K13993) and *nadM* (K00952) functions, with paralogous families
combined at the KO level. Boxes show the median and interquartile range,
whiskers extend to 1.5 times the interquartile range, and points represent
individual genomes. HighT, n = 11; LowT, n = 23; other *Nitrosopumilus*, n = 74.
"""
    (OUT / "FIGURE_LEGEND_DRAFT.md").write_text(text, encoding="utf-8")


def audit_rgb_exports(pdf, svg, png, audit_colors, expected_cell_colors):
    svg_text = svg.read_text(encoding="utf-8").lower()
    if "opacity" in svg_text:
        raise RuntimeError("SVG contains an opacity attribute.")
    if "<text" not in svg_text:
        raise RuntimeError("SVG text was not preserved as editable text.")
    if svg_text.count("<image") > 1:
        raise RuntimeError(
            "SVG contains an unexpected raster image beyond the colorbar."
        )
    for hex_color, expected_count in expected_cell_colors.items():
        if svg_text.count(hex_color.lower()) < expected_count:
            raise RuntimeError(
                f"SVG does not contain {expected_count} independent cells "
                f"with color {hex_color}."
            )

    image = Image.open(png)
    rgb = image.convert("RGB")
    color_counts = rgb.getcolors(maxcolors=rgb.width * rgb.height)
    color_counts = {color: count for count, color in color_counts}
    if "A" in image.mode and image.getchannel("A").getextrema() != (255, 255):
        raise RuntimeError("PNG contains non-opaque pixels.")

    pdf_raw = pdf.read_bytes()
    if b"/DeviceCMYK" in pdf_raw:
        raise RuntimeError("PDF contains a DeviceCMYK color space.")
    if b"/DeviceRGB" not in pdf_raw:
        raise RuntimeError("PDF does not declare a DeviceRGB color space.")
    decoded_streams = []
    for match in re.finditer(
        rb"stream\r?\n(.*?)\r?\nendstream", pdf_raw, flags=re.DOTALL
    ):
        try:
            decoded_streams.append(zlib.decompress(match.group(1)))
        except zlib.error:
            continue
    decoded = b"\n".join(decoded_streams)
    pdf_rgb_values = [
        tuple(float(value) for value in match)
        for match in re.findall(
            rb"([0-9]*\.?[0-9]+)\s+([0-9]*\.?[0-9]+)\s+"
            rb"([0-9]*\.?[0-9]+)\s+(?:rg|RG)",
            decoded,
        )
    ]
    pdf_gray_values = [
        float(value)
        for value in re.findall(
            rb"(?<![0-9.])([0-9]*\.?[0-9]+)\s+(?:g|G)(?![A-Za-z])",
            decoded,
        )
    ]

    for name, hex_color in audit_colors.items():
        normalized = hex_color.lower()
        if normalized not in svg_text:
            raise RuntimeError(f"SVG is missing {name}: {hex_color}")
        rgb_255 = tuple(
            int(normalized[index:index + 2], 16) for index in (1, 3, 5)
        )
        if color_counts.get(rgb_255, 0) == 0:
            raise RuntimeError(f"PNG is missing {name}: {hex_color}")
        rgb_unit = tuple(value / 255 for value in rgb_255)
        rgb_match = any(
            max(abs(observed - expected) for observed, expected in zip(
                pdf_value, rgb_unit
            )) < 0.002
            for pdf_value in pdf_rgb_values
        )
        gray_match = (
            rgb_255[0] == rgb_255[1] == rgb_255[2]
            and any(abs(value - rgb_unit[0]) < 0.002 for value in pdf_gray_values)
        )
        if not (rgb_match or gray_match):
            closest = min(
                pdf_rgb_values,
                key=lambda pdf_value: max(
                    abs(observed - expected)
                    for observed, expected in zip(pdf_value, rgb_unit)
                ),
                default=None,
            )
            raise RuntimeError(
                f"PDF is missing {name}: {hex_color}; closest RGB={closest}"
            )

    print(
        "RGB audit passed: SVG has editable text, independent vector cells, "
        "and no opacity; "
        "PNG is fully opaque; PDF uses DeviceRGB with no DeviceCMYK."
    )
    print("Audited colors: " + ", ".join(dict.fromkeys(audit_colors.values())))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    family = load_matrix(FAMILY_MATRIX)
    ko = load_matrix(KO_MATRIX)
    source_a = prevalence_data(family, PANEL_A, "a")
    source_b = prevalence_data(family, PANEL_B, "b")
    source_c = copy_number_data(ko)
    source = pd.concat([source_a, source_b, source_c], ignore_index=True)
    source.to_csv(
        OUT / "Fig5_source_data.tsv", sep="\t", index=False,
        float_format="%.6g",
    )

    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 7.5,
        "axes.linewidth": 0.8,
        "figure.facecolor": BACKGROUND_COLOR,
        "axes.facecolor": BACKGROUND_COLOR,
        "savefig.facecolor": BACKGROUND_COLOR,
        "text.color": TEXT_COLOR,
        "axes.edgecolor": TEXT_COLOR,
        "axes.labelcolor": TEXT_COLOR,
        "xtick.color": TEXT_COLOR,
        "ytick.color": TEXT_COLOR,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.transparent": False,
    })
    fig = plt.figure(figsize=(7.2, 3.08), constrained_layout=False)
    grid = fig.add_gridspec(
        1, 3, width_ratios=[4.60, 0.72, 3.00],
        left=0.075, right=0.99, top=0.78, bottom=0.30, wspace=0.27,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1], sharey=ax_a)
    ax_c = fig.add_subplot(grid[0, 2])

    image = draw_prevalence(
        ax_a, source_a, PANEL_A, "LowT depletion", "a",
        functional_groups=PANEL_A_GROUPS,
    )
    draw_prevalence(
        ax_b, source_b, PANEL_B, "LowT-biased\ncandidate", "b",
        highlight=(1, 0),
    )
    ax_b.tick_params(labelleft=False)
    draw_copy_number(ax_c, source_c)

    colorbar_ax = fig.add_axes([0.18, 0.075, 0.31, 0.027])
    colorbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=Normalize(0, 100), cmap=image.cmap),
        cax=colorbar_ax, orientation="horizontal",
    )
    colorbar.set_ticks([0, 50, 100])
    colorbar.ax.tick_params(labelsize=6.5, length=2, pad=1)
    colorbar.set_label("Prevalence (%)", fontsize=7, labelpad=2)
    colorbar.outline.set_linewidth(0.6)

    output_base = OUT / "Fig5_asymmetric_accessory_repertoires_RGB"
    pdf = output_base.with_suffix(".pdf")
    svg = output_base.with_suffix(".svg")
    png = output_base.with_suffix(".png")
    save_options = {
        "facecolor": BACKGROUND_COLOR,
        "edgecolor": BACKGROUND_COLOR,
        "transparent": False,
    }
    fig.savefig(pdf, **save_options)
    fig.savefig(svg, **save_options)
    fig.savefig(png, dpi=600, **save_options)
    plt.close(fig)
    write_legend()
    prevalence_values = pd.concat(
        [source_a["value"], source_b["value"]], ignore_index=True
    )
    expected_cell_colors = Counter(
        prevalence_hex(value) for value in prevalence_values
    )
    cell_audit_colors = {
        f"prevalence_cell_{index}": hex_color
        for index, hex_color in enumerate(expected_cell_colors)
    }
    audit_rgb_exports(
        pdf, svg, png,
        {**AUDIT_COLORS, **cell_audit_colors},
        expected_cell_colors,
    )

    summary = (
        source_c.groupby(["feature_label", "group"])["value"]
        .agg(["count", "mean", "median", "min", "max"])
    )
    print(summary.to_string())
    print(f"PDF: {pdf}")
    print(f"SVG: {svg}")
    print(f"PNG: {png}")


if __name__ == "__main__":
    main()
