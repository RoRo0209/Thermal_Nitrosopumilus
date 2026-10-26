#!/usr/bin/env python
"""Reanalyse inStrain metrics at the sample-by-ecotype level.

This script deliberately excludes the legacy ``pN_pS`` column because it was
calculated as a nonsynonymous/synonymous SNV count ratio without correcting for
the numbers of nonsynonymous and synonymous sites.
"""

# Release path adapter: data and outputs live in a separate user workspace.
import os as _release_os
from pathlib import Path as _ReleasePath
if not _release_os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the analysis workspace")
_release_root = _ReleasePath(_release_os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
_release_os.chdir(_release_root / 'instrain')


from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon


INPUT = Path("out/instrain_coverage_dayu1.xlsx")
OUTDIR = Path("reanalysis")
MIN_COVERAGE = 10.0
MIN_BREADTH = 0.5
MIN_LINKED_SNVS = 100


def weighted_mean(values, weights):
    values = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    weights = pd.to_numeric(weights, errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not keep.any():
        return np.nan
    return np.average(values[keep], weights=weights[keep])


def aggregate_group(group):
    length = pd.to_numeric(group["length"], errors="coerce")
    snvs = pd.to_numeric(group["SNV_count"], errors="coerce")
    linked = pd.to_numeric(group["linked_SNV_count"], errors="coerce")
    return pd.Series(
        {
            "n_genomes": len(group),
            "coverage_median": group["coverage"].median(),
            "breadth_median": group["breadth"].median(),
            "SNVs_per_kbp": snvs.sum() / (length.sum() / 1000),
            "nucleotide_diversity": weighted_mean(
                group["nucl_diversity_rarefied"], length
            ),
            "D_prime": weighted_mean(group["d_prime_mean"], linked),
            "linked_SNV_count": linked.sum(),
        }
    )


def main():
    OUTDIR.mkdir(exist_ok=True)
    df = pd.read_excel(INPUT)
    df["clade"] = df["clade"].str.strip().str.lower()

    numeric = [
        "coverage",
        "breadth",
        "length",
        "SNV_count",
        "nucl_diversity_rarefied",
        "d_prime_mean",
        "linked_SNV_count",
        "pN_pS",
    ]
    for column in numeric:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    qc = df.loc[
        df["clade"].isin(["hight", "lowt"])
        & (df["coverage"] >= MIN_COVERAGE)
        & (df["breadth"] >= MIN_BREADTH)
    ].copy()

    sample_level = (
        qc.groupby(["sample", "clade"], observed=True)
        .apply(aggregate_group, include_groups=False)
        .reset_index()
    )
    sample_level.loc[
        sample_level["linked_SNV_count"] < MIN_LINKED_SNVS, "D_prime"
    ] = np.nan
    sample_level.to_csv(
        OUTDIR / "sample_ecotype_metrics_qc.tsv", sep="\t", index=False
    )

    records = []
    for metric in ["SNVs_per_kbp", "nucleotide_diversity", "D_prime"]:
        high = sample_level.loc[sample_level.clade == "hight", metric].dropna()
        low = sample_level.loc[sample_level.clade == "lowt", metric].dropna()
        if len(high) and len(low):
            _, p_unpaired = mannwhitneyu(high, low, alternative="two-sided")
        else:
            p_unpaired = np.nan

        paired = sample_level.pivot(index="sample", columns="clade", values=metric)
        paired = paired.dropna(subset=["hight", "lowt"])
        if len(paired) >= 3 and (paired.hight != paired.lowt).any():
            _, p_paired = wilcoxon(paired.hight, paired.lowt)
        else:
            p_paired = np.nan

        rho, p_cov = spearmanr(
            sample_level[metric],
            np.log10(sample_level["coverage_median"]),
            nan_policy="omit",
        )
        records.append(
            {
                "metric": metric,
                "HighT_n": len(high),
                "LowT_n": len(low),
                "HighT_median": high.median(),
                "LowT_median": low.median(),
                "Mann_Whitney_p": p_unpaired,
                "paired_samples_n": len(paired),
                "paired_Wilcoxon_p": p_paired,
                "coverage_Spearman_rho": rho,
                "coverage_Spearman_p": p_cov,
            }
        )

    stats = pd.DataFrame(records)
    stats.to_csv(OUTDIR / "sample_ecotype_metric_tests.tsv", sep="\t", index=False)

    legacy = df["pN_pS"].dropna()
    audit = [
        "Population-genomic metric reanalysis",
        "====================================",
        f"Input: {INPUT}",
        f"QC: coverage >= {MIN_COVERAGE}x; breadth >= {MIN_BREADTH}",
        "Statistical unit: sample x ecotype after within-sample aggregation",
        f"Rows before QC: {len(df)}",
        f"Rows after QC: {len(qc)}",
        f"Sample-ecotype observations: {len(sample_level)}",
        "",
        "Important pN/pS audit:",
        "The legacy pN_pS column was not analysed. It was generated from",
        "(N_SNV_count + 0.5) / (S_SNV_count + 0.5), without normalization by",
        "N_sites and S_sites. Standard pN/pS must be recalculated from gene_info",
        "as (sum(N_SNVs)/sum(N_sites)) / (sum(S_SNVs)/sum(S_sites)).",
        f"Legacy non-null values retained only for traceability: {len(legacy)}",
        "",
        stats.to_string(index=False),
    ]
    (OUTDIR / "population_metric_audit.txt").write_text(
        "\n".join(audit), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
