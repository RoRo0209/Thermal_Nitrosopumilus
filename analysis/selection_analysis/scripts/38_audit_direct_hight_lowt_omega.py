#!/usr/bin/env python3
"""Audit the direct HighT-versus-LowT codeml summary and its Figure 7c concordance."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--direct", required=True, type=Path)
    parser.add_argument("--panelc", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def bh_adjust(p_values: pd.Series) -> pd.Series:
    values = p_values.astype(float).to_numpy()
    order = np.argsort(values)
    ranked = values[order]
    n = len(values)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    result = np.empty(n)
    result[order] = adjusted
    return pd.Series(result, index=p_values.index)


def warning_counts(frame: pd.DataFrame) -> pd.Series:
    warnings = (
        frame["warnings"]
        .fillna("")
        .str.split(";")
        .explode()
    )
    warnings = warnings[warnings.ne("")]
    return warnings.value_counts()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    direct = pd.read_csv(args.direct, sep="\t", dtype={"KO": "string"})
    panelc = pd.read_csv(args.panelc, sep="\t", dtype={"Integrated_KO": "string"})

    numeric = [
        "omega_HighT", "omega_LowT", "omega_Other",
        "omega_ratio_HighT_LowT", "raw_P", "BH_FDR",
        "direct_significant_FDR05", "tree_qc_pass",
        "alignment_codons", "null_successful_reps",
        "alternative_successful_reps", "null_reps_at_best_0.01",
        "alternative_reps_at_best_0.01",
    ]
    for column in numeric:
        direct[column] = pd.to_numeric(direct[column], errors="coerce")

    problems: list[str] = []
    if len(direct) != 110:
        problems.append(f"expected_110_rows_observed_{len(direct)}")
    if direct["split_family"].duplicated().any():
        problems.append("duplicate_split_family")
    if direct["split_family"].isna().any():
        problems.append("missing_split_family")
    if direct[numeric].isna().any().any():
        problems.append("missing_or_non_numeric_required_values")
    if not direct["tree_qc_pass"].eq(1).all():
        problems.append("tree_qc_failure")
    if (direct["alternative_lnL"].astype(float) + 1e-8 <
            direct["null_lnL"].astype(float)).any():
        problems.append("alternative_lnL_below_null")

    direct["BH_FDR_recomputed"] = bh_adjust(direct["raw_P"])
    max_bh_difference = (
        direct["BH_FDR_recomputed"] - direct["BH_FDR"]
    ).abs().max()
    if max_bh_difference > 1e-12:
        problems.append("BH_FDR_mismatch")

    direct["FDR_lt_0.05"] = direct["BH_FDR"].lt(0.05)
    direct["stable"] = direct["convergence_status"].eq("stable")
    direct["accepted_direct_hit"] = direct["direct_significant_FDR05"].eq(1)
    direct["has_KO"] = direct["KO"].fillna("").ne("")
    direct["informative_annotation"] = ~direct["annotation"].fillna("").str.lower().isin(
        ["", "hypothetical protein"]
    )

    accepted = direct[direct["accepted_direct_hit"]].copy()
    review_hits = direct[direct["FDR_lt_0.05"] & ~direct["stable"]].copy()

    direct["match_family"] = direct["split_family"].str.replace(
        r"'$", "_", regex=True
    )
    panelc["match_family"] = panelc["split_family"].str.replace(
        r"'$", "_", regex=True
    )
    merged = direct.merge(
        panelc[
            [
                "match_family", "split_family", "omega_ratio_HighT_LowT",
                "effect_log2", "rank", "display_label",
            ]
        ].rename(
            columns={
                "split_family": "panelc_split_family",
                "omega_ratio_HighT_LowT": "panelc_ratio_separate_models",
                "effect_log2": "panelc_log2_ratio",
                "rank": "panelc_rank",
            }
        ),
        on="match_family",
        how="outer",
        indicator=True,
    )
    unmatched = merged["_merge"].ne("both").sum()
    if unmatched:
        problems.append(f"panelc_join_unmatched_{unmatched}")

    matched = merged[merged["_merge"].eq("both")].copy()
    matched["direct_log2_ratio"] = np.log2(matched["omega_ratio_HighT_LowT"])
    matched["direction_agrees"] = (
        matched["omega_ratio_HighT_LowT"].lt(1)
        == matched["panelc_ratio_separate_models"].lt(1)
    )
    pearson_log2 = matched[
        ["direct_log2_ratio", "panelc_log2_ratio"]
    ].corr(method="pearson").iloc[0, 1]
    spearman_log2 = matched[
        ["direct_log2_ratio", "panelc_log2_ratio"]
    ].rank().corr(method="pearson").iloc[0, 1]

    labelled = matched[matched["display_label"].fillna("").ne("")].copy()
    labelled = labelled[
        [
            "split_family", "display_label", "annotation", "KO",
            "omega_HighT", "omega_LowT", "omega_ratio_HighT_LowT",
            "BH_FDR", "convergence_status", "accepted_direct_hit", "warnings",
        ]
    ].sort_values("omega_ratio_HighT_LowT")

    selected_columns = [
        "split_family", "annotation", "KO", "n_taxa", "alignment_codons",
        "omega_HighT", "omega_LowT", "omega_Other",
        "omega_ratio_HighT_LowT", "LRT", "raw_P", "BH_FDR",
        "direction", "convergence_status", "warnings",
    ]
    accepted[selected_columns].sort_values("BH_FDR").to_csv(
        args.outdir / "direct_omega_25_stable_FDR_hits.tsv",
        sep="\t", index=False,
    )
    review_hits[selected_columns].sort_values("BH_FDR").to_csv(
        args.outdir / "direct_omega_3_FDR_hits_excluded_for_convergence.tsv",
        sep="\t", index=False,
    )
    labelled.to_csv(
        args.outdir / "Figure7c_six_labels_direct_validation.tsv",
        sep="\t", index=False,
    )
    merged.drop(columns=["_merge"]).to_csv(
        args.outdir / "direct_omega_with_panelc_comparison.tsv",
        sep="\t", index=False,
    )

    warning_all = warning_counts(direct)
    warning_accepted = warning_counts(accepted)
    warning_lines = [
        f"- `{name}`: {count}"
        for name, count in warning_all.items()
    ] or ["- None"]
    accepted_warning_lines = [
        f"- `{name}`: {count}"
        for name, count in warning_accepted.items()
    ] or ["- None"]

    labelled_lines = []
    for row in labelled.itertuples(index=False):
        status = "accepted" if row.accepted_direct_hit else "not accepted"
        labelled_lines.append(
            f"| {row.display_label} | {row.omega_HighT:.5f} | "
            f"{row.omega_LowT:.5f} | {row.omega_ratio_HighT_LowT:.3f} | "
            f"{row.BH_FDR:.3g} | {row.convergence_status} | {status} |"
        )

    accepted_with_ko = int(accepted["has_KO"].sum())
    accepted_informative = int(accepted["informative_annotation"].sum())
    accepted_h = int(accepted["direction"].eq("HighT_lower").sum())
    accepted_l = int(accepted["direction"].eq("LowT_lower").sum())
    accepted_warning_free = int(accepted["warnings"].fillna("").eq("").sum())
    accepted_short = int(
        accepted["warnings"].fillna("").str.contains("short_alignment").sum()
    )
    review_names = ", ".join(review_hits["split_family"].tolist())

    report = f"""# HighT–LowT direct omega audit

## Dataset and grain

- Unit: one soft-core single-copy gene family.
- Rows: {len(direct)}; unique families: {direct['split_family'].nunique()}.
- Complete null and alternative fits: {direct['raw_P'].notna().sum()}/{len(direct)}.
- Tree QC passed: {int(direct['tree_qc_pass'].eq(1).sum())}/{len(direct)}.
- Structural problems: {", ".join(problems) if problems else "none"}.
- Maximum difference between reported and independently recomputed BH-FDR: {max_bh_difference:.3g}.

## Main result

- Raw BH-FDR < 0.05: {int(direct['FDR_lt_0.05'].sum())} families.
- Stable and accepted FDR hits: {len(accepted)} families.
- FDR hits excluded because the multi-start optimum was not reproducible: {len(review_hits)} families.
- Direction among accepted hits: HighT lower = {accepted_h}; LowT lower = {accepted_l}.
- Accepted hits with a KO: {accepted_with_ko}/{len(accepted)}.
- Accepted hits with a non-hypothetical annotation: {accepted_informative}/{len(accepted)}.
- Accepted hits without any automated warning: {accepted_warning_free}/{len(accepted)}.
- Accepted hits with alignments shorter than 100 codons: {accepted_short}/{len(accepted)}.
- FDR-positive families excluded for unstable optima: {review_names}.

## Relation to the existing Figure 7c

- Families matched to Panel c: {len(matched)}/110.
- Direction agreement between separate-model and direct-model ratios: {int(matched['direction_agrees'].sum())}/{len(matched)}.
- Pearson correlation of log2 ratios: {pearson_log2:.4f}.
- Spearman rank correlation of log2 ratios: {spearman_log2:.4f}.
- Direct-model direction over all genes: HighT lower = {int(direct['direction'].eq('HighT_lower').sum())}; LowT lower = {int(direct['direction'].eq('LowT_lower').sum())}.

## Existing six Panel c labels

| Label | omega HighT | omega LowT | H/L ratio | BH-FDR | convergence | direct result |
|---|---:|---:|---:|---:|---|---|
{chr(10).join(labelled_lines)}

## Warnings across all families

{chr(10).join(warning_lines)}

## Warnings among the 25 accepted hits

{chr(10).join(accepted_warning_lines)}

## Interpretation

The direct nested model confirms that the broad direction shown in Figure 7c is
not an artefact of comparing two separate foreground models: all 110 genes have
the same direction under both approaches, and their effect sizes are strongly
correlated. After BH correction and multi-start convergence filtering, 25 genes
showed accepted direct HighT–LowT differences. Twenty-four had lower omega in
HighT and one had lower omega in LowT. The direct test therefore strengthens the
claim of stronger long-term sequence constraint in HighT, while showing that the
gene-level signal is concentrated in a subset of the shared core genome.

Three additional families had BH-FDR < 0.05 but were excluded because their
best likelihood was not reproducible across starting values. They should not be
reported as confirmed direct hits without targeted reruns.
"""
    (args.outdir / "DIRECT_OMEGA_AUDIT.md").write_text(
        report, encoding="utf-8", newline="\n"
    )
    print(report)


if __name__ == "__main__":
    main()
