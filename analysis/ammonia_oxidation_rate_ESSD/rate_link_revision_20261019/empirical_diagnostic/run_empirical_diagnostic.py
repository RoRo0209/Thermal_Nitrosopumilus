"""Diagnostic only: matched modeled relative abundance versus observed AO rate.

No figures, hyperparameter search, future extrapolation, or input modifications.
Run with Python containing numpy and pandas; scipy is not required.
"""
from pathlib import Path
import hashlib
import json
import platform
import sys

import numpy as np
import pandas as pd


import os
if not os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the separate analysis workspace")
_workspace = Path(os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
OUT = _workspace / 'ammonia_oxidation_rate_ESSD/rate_link_revision_20261019/empirical_diagnostic'
OUT.mkdir(parents=True, exist_ok=True)
AO = OUT.parents[1]
SOURCE = AO / "fixed_coefficient_scenario_20261016"
INPUTS = [SOURCE / "source_calibration_39_grids.tsv",
          SOURCE / "source_included_100_records.tsv"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(frame, filename):
    frame.to_csv(OUT / filename, sep="\t", index=False, float_format="%.15g")


def fit(x, y):
    design = np.column_stack([np.ones(len(x)), x])
    assert np.linalg.matrix_rank(design) == 2
    return np.linalg.lstsq(design, y, rcond=None)[0]


def main():
    hashes = {str(p): digest(p) for p in INPUTS}
    raw_pairs = pd.read_csv(INPUTS[0], sep="\t")
    records = pd.read_csv(INPUTS[1], sep="\t")
    assert len(raw_pairs) == 39 and len(records) == 100
    assert raw_pairs.grid_id.is_unique
    assert not raw_pairs[["Longitude", "Latitude"]].duplicated().any()
    assert records.excel_row.is_unique
    assert records.study.notna().all() and records.study.nunique() == 18
    assert set(records.grid_id) == set(raw_pairs.grid_id)
    assert records.rate.gt(0).all()
    study_counts = records.groupby("grid_id").study.nunique()
    assert study_counts.eq(1).all(), "Each grid must belong to one study"
    assert raw_pairs.n_studies.eq(1).all()
    rebuilt = records.groupby("grid_id", as_index=False).agg(
        study=("study", "first"), reconstructed_rate=("rate", "median"),
        reconstructed_n_records=("rate", "size"))
    pairs = raw_pairs[["grid_id", "Longitude", "Latitude", "rate_median",
                       "n_records", "low_current", "high_current", "S_current"]].merge(
                           rebuilt, on="grid_id", validate="one_to_one")
    assert np.allclose(pairs.rate_median, pairs.reconstructed_rate, rtol=1e-12, atol=1e-12)
    assert pairs.n_records.eq(pairs.reconstructed_n_records).all()
    assert np.allclose(pairs.S_current, pairs.low_current + pairs.high_current,
                       rtol=1e-12, atol=1e-14)
    assert np.isfinite(pairs[["Longitude", "Latitude", "rate_median", "S_current"]]).all().all()
    assert pairs.rate_median.gt(0).all() and pairs.S_current.gt(0).all()
    pairs = pairs.rename(columns={"rate_median": "observed_grid_median_rate",
                                  "low_current": "modeled_low_current_relative_abundance",
                                  "high_current": "modeled_high_current_relative_abundance",
                                  "S_current": "modeled_total_current_relative_abundance"})
    pairs["block_lon_lower"] = np.floor(pairs.Longitude / 5).astype(int) * 5
    pairs["block_lat_lower"] = np.floor(pairs.Latitude / 5).astype(int) * 5
    pairs["spatial_5deg_block"] = pairs.block_lon_lower.astype(str) + ":" + pairs.block_lat_lower.astype(str)
    pairs["log10_modeled_total_current_relative_abundance"] = np.log10(pairs.modeled_total_current_relative_abundance)
    pairs["log10_observed_grid_median_rate"] = np.log10(pairs.observed_grid_median_rate)
    pairs["abundance_source_sha256"] = hashes[str(INPUTS[0])]
    pairs["records_source_sha256"] = hashes[str(INPUTS[1])]
    save(pairs, "cell_level_pairs.tsv")
    trace = records.merge(pairs[["grid_id", "Longitude", "Latitude", "spatial_5deg_block"]],
                          on="grid_id", validate="many_to_one")
    trace["source_file"] = str(INPUTS[1])
    trace["source_sha256"] = hashes[str(INPUTS[1])]
    save(trace, "record_to_cell_source_trace.tsv")

    x = pairs.log10_modeled_total_current_relative_abundance.to_numpy()
    y = pairs.log10_observed_grid_median_rate.to_numpy()
    s = pairs.modeled_total_current_relative_abundance.to_numpy()
    r = pairs.observed_grid_median_rate.to_numpy()
    coefficients, heldout, fold_audit = [], [], []
    full_beta = fit(x, y)
    full_pred = full_beta[0] + full_beta[1] * x
    coefficients.append({"validation": "full_sample_descriptive_only", "heldout_group": "none",
        "model": "loglog_ols", "n_train_cells": len(pairs), "n_test_cells": 0,
        "intercept_log10": full_beta[0], "slope": full_beta[1],
        "k_original_scale": np.nan})
    coefficients.append({"validation": "full_sample_descriptive_only", "heldout_group": "none",
        "model": "proportional_sum_ratio", "n_train_cells": len(pairs), "n_test_cells": 0,
        "intercept_log10": np.log10(r.sum() / s.sum()), "slope": 1.0,
        "k_original_scale": r.sum() / s.sum()})
    for validation, column in [("leave_one_study_out", "study"),
                               ("leave_one_5deg_block_out", "spatial_5deg_block")]:
        for group in sorted(pairs[column].unique()):
            test = (pairs[column] == group).to_numpy()
            train = ~test
            assert not set(pairs.loc[train, "grid_id"]) & set(pairs.loc[test, "grid_id"])
            beta = fit(x[train], y[train])
            k = r[train].sum() / s[train].sum()
            predictions = {
                "loglog_ols": beta[0] + beta[1] * x[test],
                "fold_training_logmean": np.full(test.sum(), y[train].mean()),
                "proportional_sum_ratio": np.log10(k * s[test]),
            }
            for name, intercept, slope, kval in [
                    ("loglog_ols", beta[0], beta[1], np.nan),
                    ("fold_training_logmean", y[train].mean(), 0.0, np.nan),
                    ("proportional_sum_ratio", np.log10(k), 1.0, k)]:
                coefficients.append({"validation": validation, "heldout_group": group,
                    "model": name, "n_train_cells": int(train.sum()), "n_test_cells": int(test.sum()),
                    "intercept_log10": intercept, "slope": slope, "k_original_scale": kval})
            fold_audit.append({"validation": validation, "heldout_group": group,
                "n_train_cells": int(train.sum()), "n_test_cells": int(test.sum()),
                "train_test_grid_overlap": 0,
                "test_grid_ids": ";".join(pairs.loc[test, "grid_id"].astype(str)),
                "test_log_abundance_outside_training_range_count": int(((x[test] < x[train].min()) |
                                                                       (x[test] > x[train].max())).sum())})
            for name, pred in predictions.items():
                for index, prediction in zip(np.where(test)[0], pred):
                    heldout.append({"validation": validation, "heldout_group": group,
                        "model": name, "grid_id": int(pairs.iloc[index].grid_id),
                        "study": pairs.iloc[index].study,
                        "spatial_5deg_block": pairs.iloc[index].spatial_5deg_block,
                        "observed_log10_rate": y[index], "predicted_log10_rate": prediction,
                        "observed_rate": r[index], "predicted_rate_geometric_scale": 10 ** prediction,
                        "fold_training_logmean": y[train].mean(),
                        "n_train_cells": int(train.sum()), "n_test_cells": int(test.sum())})
    heldout = pd.DataFrame(heldout)
    metrics = []
    for (validation, model), part in heldout.groupby(["validation", "model"], sort=False):
        assert len(part) == 39 and part.grid_id.is_unique
        residual = part.predicted_log10_rate - part.observed_log10_rate
        baseline_residual = part.fold_training_logmean - part.observed_log10_rate
        sse = float(np.sum(residual ** 2))
        baseline_sse = float(np.sum(baseline_residual ** 2))
        sst = float(np.sum((part.observed_log10_rate - part.observed_log10_rate.mean()) ** 2))
        fold_mse = part.assign(squared_error=residual ** 2).groupby("heldout_group").squared_error.mean()
        metrics.append({"validation": validation, "model": model, "n_cells": len(part),
            "n_groups": part.heldout_group.nunique(), "rmse_log10": np.sqrt(sse / len(part)),
            "mae_log10": float(np.mean(np.abs(residual))),
            "pooled_R2_log10": 1 - sse / sst,
            "skill_vs_fold_training_logmean": 1 - sse / baseline_sse,
            "equal_group_rmse_log10": np.sqrt(fold_mse.mean()),
            "median_multiplicative_error": 10 ** float(np.median(np.abs(residual))),
            "fraction_within_factor_2": float(np.mean(np.abs(residual) <= np.log10(2)))})
    save(pd.DataFrame(coefficients), "coefficients_summary.tsv")
    save(heldout, "heldout_predictions.tsv")
    save(pd.DataFrame(fold_audit), "fold_audit.tsv")
    save(pd.DataFrame(metrics), "validation_metrics.tsv")
    summary = {
        "status": "empirical_diagnostic_only_no_future_projection",
        "analysis_date": "2026-10-19", "input_sha256": hashes,
        "script_sha256": digest(Path(__file__)),
        "runtime": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
                    "platform": platform.platform()},
        "source_checks": {"n_cells": len(pairs), "n_records": len(records),
            "n_studies": int(pairs.study.nunique()), "n_5deg_blocks": int(pairs.spatial_5deg_block.nunique()),
            "all_cells_exactly_one_study": True, "grid_rate_medians_reproduced": True,
            "modeled_total_abundance_equals_low_plus_high": True,
            "no_cell_overlap_in_any_train_test_split": True},
        "prespecified_model": "log10(rate) = intercept + slope * log10(modeled LowT+HighT current relative abundance)",
        "fitting": "Unweighted ordinary least squares using numpy.linalg.lstsq; free intercept and slope; no tuning.",
        "comparators": ["Fold-training mean log10(rate)",
                        "k=sum(training rates)/sum(training abundances); prediction=k*heldout abundance"],
        "full_sample_descriptive_fit": {"intercept_log10": float(full_beta[0]),
            "slope": float(full_beta[1]),
            "in_sample_R2_log10": float(1 - np.sum((full_pred-y)**2) / np.sum((y-y.mean())**2)),
            "pearson_log10": float(np.corrcoef(x, y)[0, 1]),
            "spearman": float(pd.Series(x).rank().corr(pd.Series(y).rank()))},
        "metrics": metrics,
        "limitations": [
            "Abundance is modeled climate-average relative abundance, not observed contemporaneous absolute AOA abundance.",
            "Spatially matched rate observations and abundance predictions may mismatch in time, measurement, ecology, and environmental statistics.",
            "Environmental ammonia oxidation is not established as exclusively attributable to the LowT/HighT lineages.",
            "Only 39 cells from 18 studies; all records in a cell originate from one study; within-cell median pools sampling times.",
            "Validation holds out complete studies or 5-degree blocks but cannot remove upstream abundance-model training leakage if present; upstream model training membership was not evaluated here.",
            "These tests assess this specific matched-data empirical link, not whether a contemporaneous absolute-abundance/rate relationship exists in nature.",
            "Spatial block origin is fixed at multiples of 5 degrees; no blocked-design or hyperparameter search was performed.",
            "10**predicted log10 rate is on a geometric scale; no retransformation correction to arithmetic mean rate was made.",
            "No uncertainty interval or significance claim is inferred from the descriptive full-sample slope.",
            "No new figures, future projections, or changes to existing analyses were made."],
    }
    after = {str(p): digest(p) for p in INPUTS}
    assert after == hashes, "Input changed during run"
    summary["input_hashes_unchanged_after_run"] = True
    (OUT / "diagnostic_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({"full_sample_descriptive_fit": summary["full_sample_descriptive_fit"],
                      "source_checks": summary["source_checks"], "metrics": metrics}, indent=2))


if __name__ == "__main__":
    main()
