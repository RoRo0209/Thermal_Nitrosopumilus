#!/usr/bin/env python3
"""Summarize paired LowT and HighT SSP5-8.5 abundance-change directions."""

import os
from pathlib import Path
_thermal_root = os.environ.get("THERMAL_PROJECT_ROOT", "")
if not _thermal_root or not Path(_thermal_root).is_dir():
    raise SystemExit("Set THERMAL_PROJECT_ROOT to an existing, separate analysis workspace.")
THERMAL_PROJECT_ROOT = Path(_thermal_root).resolve()


import numpy as np
import pandas as pd


BASE = THERMAL_PROJECT_ROOT / 'SSP585'
TAG = "GLOBAL_nearshore_pred_delta_nearshore200km_KEEPCHL_SSP585_2080-2100.tsv"


def main() -> None:
    low = pd.read_csv(
        BASE / "NP_lowT" / "AOA_future_projection_outputs" / "GLOBAL_NEARSHORE_MAPS" / TAG,
        sep="\t",
    )[["Longitude", "Latitude", "delta_future_minus_current"]].rename(
        columns={"delta_future_minus_current": "LowT_delta_log10"}
    )
    high = pd.read_csv(
        BASE / "NP_highT" / "AOA_future_projection_outputs" / "GLOBAL_NEARSHORE_MAPS" / TAG,
        sep="\t",
    )[["Longitude", "Latitude", "delta_future_minus_current"]].rename(
        columns={"delta_future_minus_current": "HighT_delta_log10"}
    )
    data = low.merge(high, on=["Longitude", "Latitude"], how="inner")
    data = data.sort_values(["Latitude", "Longitude"]).reset_index(drop=True)
    data["latitude_sorted_index"] = np.arange(1, len(data) + 1)
    data["direction_pair"] = np.select(
        [
            (data["LowT_delta_log10"] > 0) & (data["HighT_delta_log10"] > 0),
            (data["LowT_delta_log10"] < 0) & (data["HighT_delta_log10"] < 0),
            (data["LowT_delta_log10"] > 0) & (data["HighT_delta_log10"] < 0),
        ],
        ["both_increase", "both_decrease", "LowT_up_HighT_down"],
        default="LowT_down_HighT_up",
    )

    threshold = 0.2
    for ecotype in ["LowT", "HighT"]:
        x = data[f"{ecotype}_delta_log10"]
        data[f"{ecotype}_threshold_class"] = np.select(
            [x > threshold, x < -threshold],
            ["increase", "decrease"],
            default="stable",
        )

    counts = data["direction_pair"].value_counts().rename_axis("direction_pair").reset_index(name="grid_cells")
    counts["percent"] = counts["grid_cells"] / len(data) * 100
    table = pd.crosstab(
        data["LowT_threshold_class"],
        data["HighT_threshold_class"],
    )
    data.to_csv(BASE / "LowT_HighT_latitude_sorted_change.tsv", sep="\t", index=False)
    counts.to_csv(BASE / "LowT_HighT_direction_pair_summary.tsv", sep="\t", index=False)
    table.to_csv(BASE / "LowT_HighT_threshold_pair_table.tsv", sep="\t")

    opposite = (data["LowT_delta_log10"] * data["HighT_delta_log10"] < 0).sum()
    same = (data["LowT_delta_log10"] * data["HighT_delta_log10"] > 0).sum()
    sig_opp = (
        ((data["LowT_delta_log10"] > threshold) & (data["HighT_delta_log10"] < -threshold))
        | ((data["LowT_delta_log10"] < -threshold) & (data["HighT_delta_log10"] > threshold))
    ).sum()
    sig_same = (
        ((data["LowT_delta_log10"] > threshold) & (data["HighT_delta_log10"] > threshold))
        | ((data["LowT_delta_log10"] < -threshold) & (data["HighT_delta_log10"] < -threshold))
    ).sum()

    print(f"rows={len(data)}")
    print("\nANY-DIRECTION PAIRS")
    print(counts.to_string(index=False))
    print(f"\nopposite_any={opposite} ({opposite / len(data) * 100:.2f}%)")
    print(f"same_any={same} ({same / len(data) * 100:.2f}%)")
    print(f"pearson_correlation={data[['LowT_delta_log10', 'HighT_delta_log10']].corr().iloc[0, 1]:.4f}")
    print("\nTHRESHOLD +/-0.2")
    print(table.to_string())
    print(f"both_significant_opposite={sig_opp}")
    print(f"both_significant_same={sig_same}")


if __name__ == "__main__":
    main()
