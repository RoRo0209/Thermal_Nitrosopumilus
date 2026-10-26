import os
from pathlib import Path
_thermal_root = os.environ.get("THERMAL_PROJECT_ROOT", "")
if not _thermal_root or not Path(_thermal_root).is_dir():
    raise SystemExit("Set THERMAL_PROJECT_ROOT to an existing, separate analysis workspace.")
THERMAL_PROJECT_ROOT = Path(_thermal_root).resolve()


import pandas as pd


BASE = THERMAL_PROJECT_ROOT / "SSP585"
TAG = "nearshore200km_KEEPCHL_SSP585_2080-2100"
SST_CUT = 22.0
CHL_CUT = 1.5
DELTA_CUT = 0.2


def load_delta(group: str) -> pd.DataFrame:
    folder = "NP_lowT" if group == "LowT" else "NP_highT"
    path = (
        BASE
        / folder
        / "AOA_future_projection_outputs"
        / "GLOBAL_NEARSHORE_MAPS"
        / f"GLOBAL_nearshore_pred_delta_{TAG}.tsv"
    )
    df = pd.read_csv(path, sep="\t")
    out = df[["Longitude", "Latitude", "temp", "chl", "delta_future_minus_current"]].copy()
    out["group"] = group
    out = out.rename(columns={"temp": "future_sst", "chl": "future_chl", "delta_future_minus_current": "delta_log10"})
    out = out.dropna(subset=["future_sst", "future_chl", "delta_log10"])
    return out


def classify_change(delta: float) -> str:
    if delta > DELTA_CUT:
        return "Increase"
    if delta < -DELTA_CUT:
        return "Decrease"
    return "Stable"


def classify_sector(sst: float, chl: float) -> str:
    if sst >= SST_CUT and chl >= CHL_CUT:
        return "Warm + Productive"
    if sst >= SST_CUT and chl < CHL_CUT:
        return "Warm + Low-productivity"
    if sst < SST_CUT and chl >= CHL_CUT:
        return "Cool + Productive"
    return "Cool + Low-productivity"


def summarize(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = df.copy()
    df["change_class"] = df["delta_log10"].map(classify_change)
    df["environmental_sector"] = [
        classify_sector(sst, chl) for sst, chl in zip(df["future_sst"], df["future_chl"])
    ]
    order_change = ["Increase", "Stable", "Decrease"]
    order_sector = [
        "Warm + Productive",
        "Warm + Low-productivity",
        "Cool + Productive",
        "Cool + Low-productivity",
    ]
    counts = (
        df.groupby(["group", "change_class", "environmental_sector"], dropna=False)
        .size()
        .rename("n")
        .reset_index()
    )
    all_index = pd.MultiIndex.from_product(
        [df["group"].unique(), order_change, order_sector],
        names=["group", "change_class", "environmental_sector"],
    )
    counts = (
        counts.set_index(["group", "change_class", "environmental_sector"])
        .reindex(all_index, fill_value=0)
        .reset_index()
    )
    counts["total_in_class"] = counts.groupby(["group", "change_class"])["n"].transform("sum")
    counts["fraction"] = counts["n"] / counts["total_in_class"]
    counts["percent"] = counts["fraction"] * 100

    class_summary = (
        df.groupby(["group", "change_class"])
        .agg(
            n=("delta_log10", "size"),
            median_delta=("delta_log10", "median"),
            median_future_sst=("future_sst", "median"),
            median_future_chl=("future_chl", "median"),
        )
        .reset_index()
    )
    class_summary["percent_of_group"] = (
        class_summary["n"] / class_summary.groupby("group")["n"].transform("sum") * 100
    )
    return counts, class_summary


def main() -> None:
    combined = pd.concat([load_delta("LowT"), load_delta("HighT")], ignore_index=True)
    counts, class_summary = summarize(combined)

    counts.to_csv(BASE / "LowT_HighT_environmental_sector_by_change_class.tsv", sep="\t", index=False)
    class_summary.to_csv(BASE / "LowT_HighT_change_class_environment_summary.tsv", sep="\t", index=False)

    low = counts[counts["group"] == "LowT"]
    high = counts[counts["group"] == "HighT"]

    def pct(table: pd.DataFrame, change: str, sector: str) -> float:
        row = table[(table["change_class"] == change) & (table["environmental_sector"] == sector)]
        return float(row["percent"].iloc[0])

    low_inc_prod = pct(low, "Increase", "Warm + Productive") + pct(low, "Increase", "Cool + Productive")
    low_dec_warm_prod = pct(low, "Decrease", "Warm + Productive")
    low_dec_warm_low = pct(low, "Decrease", "Warm + Low-productivity")
    low_dec_warm = low_dec_warm_prod + low_dec_warm_low
    high_inc_prod = pct(high, "Increase", "Warm + Productive") + pct(high, "Increase", "Cool + Productive")
    high_dec_warm_low = pct(high, "Decrease", "Warm + Low-productivity")

    lines = [
        "# LowT environmental-context follow-up",
        "",
        f"Thresholds matched the HighT Fig. 8f-g analysis: Δlog10 abundance > {DELTA_CUT} for Increase, < -{DELTA_CUT} for Decrease, SST = {SST_CUT:g} °C and chlorophyll-a = {CHL_CUT:g} mg m-3 for environmental sectors.",
        "",
        "## Main quantitative result",
        "",
        f"- LowT Increase cells in productive sectors: {low_inc_prod:.1f}%.",
        f"- LowT Decrease cells in warm sectors: {low_dec_warm:.1f}% ({low_dec_warm_prod:.1f}% warm productive and {low_dec_warm_low:.1f}% warm low-productivity).",
        f"- HighT Increase cells in productive sectors: {high_inc_prod:.1f}%.",
        f"- HighT Decrease cells in warm low-productivity waters: {high_dec_warm_low:.1f}%.",
        "",
        "## Interpretation",
        "",
        "LowT also shows a structured environmental context, but the pattern differs from HighT. LowT increase cells are concentrated mainly in cool productive sectors, whereas LowT decrease cells are concentrated in warm sectors overall rather than specifically in warm low-productivity waters. This supports the reviewer-facing statement that LowT responses are not explained by temperature alone, while avoiding overextending the HighT-specific productivity mechanism to LowT.",
        "",
        "## Suggested manuscript wording",
        "",
        "LowT responses were also spatially heterogeneous, with regional gains and losses rather than a simple decline under warming. A parallel stratification showed that LowT increases were concentrated mainly in cool, productive waters, whereas LowT decreases were enriched across warm sectors (Supplementary Fig. Sx), indicating that temperature alone did not explain its projected trajectory.",
    ]
    (BASE / "LowT_environmental_context_followup_summary.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
