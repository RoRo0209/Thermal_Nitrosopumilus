#!/usr/bin/env python3
"""Summarize direct HighT-versus-LowT three-category codeml models."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


LNL_RE = re.compile(
    r"lnL\(.*?\):\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
)
OMEGA_RE = re.compile(r"w \(dN/dS\) for branches:\s*(.+)")
TREE_DN_RE = re.compile(r"tree length for dN:\s*([0-9.Ee+-]+)")
TREE_DS_RE = re.compile(r"tree length for dS:\s*([0-9.Ee+-]+)")
BRANCH_ROW_RE = re.compile(
    r"^\s*\S+\.\.\S+\s+"
    r"([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+"
    r"([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+"
    r"([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s*$",
    re.MULTILINE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-dir", type=Path, required=True)
    parser.add_argument("--tree-qc", type=Path, required=True)
    parser.add_argument("--annotations", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=fields, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_float_match(pattern: re.Pattern[str], text: str) -> float | None:
    match = pattern.search(text)
    return float(match.group(1)) if match else None


def parse_omegas(text: str) -> list[float]:
    match = OMEGA_RE.search(text)
    if not match:
        return []
    values: list[float] = []
    for token in match.group(1).split():
        try:
            values.append(float(token))
        except ValueError:
            continue
    return values


def parse_branch_rows(text: str) -> list[dict[str, float]]:
    rows = []
    for match in BRANCH_ROW_RE.finditer(text):
        t, n_sites, s_sites, omega, dn, ds, n_dn, s_ds = map(float, match.groups())
        rows.append({
            "t": t, "N": n_sites, "S": s_sites, "omega": omega,
            "dN": dn, "dS": ds, "N_dN": n_dn, "S_dS": s_ds,
        })
    return rows


def best_run(model_dir: Path, prefix: str) -> dict | None:
    runs = []
    for run_dir in sorted(model_dir.glob(prefix + "_rep*")):
        result = run_dir / "results.txt"
        if not result.is_file() or result.stat().st_size == 0:
            continue
        text = result.read_text(encoding="utf-8", errors="replace")
        lnl_match = LNL_RE.search(text)
        completed = "Time used:" in text
        if not lnl_match or not completed:
            continue
        lnl = float(lnl_match.group(1))
        omegas = parse_omegas(text)
        runs.append({
            "lnL": lnl,
            "run_dir": str(run_dir),
            "text": text,
            "omegas": omegas,
            "tree_dN": parse_float_match(TREE_DN_RE, text),
            "tree_dS": parse_float_match(TREE_DS_RE, text),
            "branch_rows": parse_branch_rows(text),
        })
    if not runs:
        return None
    runs.sort(key=lambda row: row["lnL"], reverse=True)
    best = runs[0]
    best["n_successful_reps"] = len(runs)
    best["lnL_rep_spread"] = runs[0]["lnL"] - runs[-1]["lnL"]
    best["n_reps_at_best_0.01"] = sum(
        runs[0]["lnL"] - row["lnL"] <= 0.01 for row in runs
    )
    best["all_rep_lnL"] = ",".join(f"{row['lnL']:.6f}" for row in runs)
    return best


def chi1_sf(statistic: float) -> float:
    return math.erfc(math.sqrt(max(statistic, 0.0) / 2.0))


def bh_adjust(values: list[float | None]) -> list[float | None]:
    indexed = [
        (index, value)
        for index, value in enumerate(values)
        if value is not None and math.isfinite(value)
    ]
    indexed.sort(key=lambda item: item[1])
    adjusted: dict[int, float] = {}
    running = 1.0
    total = len(indexed)
    for reverse_index in range(total - 1, -1, -1):
        original_index, value = indexed[reverse_index]
        rank = reverse_index + 1
        running = min(running, value * total / rank)
        adjusted[original_index] = min(running, 1.0)
    return [adjusted.get(index) for index in range(len(values))]


def class_diagnostics(
    branch_rows: list[dict[str, float]], omegas: list[float], class_index: int
) -> dict[str, float | int | None]:
    if not branch_rows or class_index >= len(omegas):
        return {
            "n_branches": 0, "sum_dN": None, "sum_dS": None,
            "median_dS": None, "zero_dS_fraction": None,
        }
    assignments = []
    for row in branch_rows:
        nearest = min(
            range(len(omegas)),
            key=lambda index: abs(row["omega"] - omegas[index]),
        )
        if nearest == class_index:
            assignments.append(row)
    if not assignments:
        return {
            "n_branches": 0, "sum_dN": 0.0, "sum_dS": 0.0,
            "median_dS": None, "zero_dS_fraction": None,
        }
    ds_values = sorted(row["dS"] for row in assignments)
    middle = len(ds_values) // 2
    median_ds = (
        ds_values[middle]
        if len(ds_values) % 2
        else (ds_values[middle - 1] + ds_values[middle]) / 2
    )
    return {
        "n_branches": len(assignments),
        "sum_dN": sum(row["dN"] for row in assignments),
        "sum_dS": sum(row["dS"] for row in assignments),
        "median_dS": median_ds,
        "zero_dS_fraction": sum(row["dS"] <= 1e-8 for row in assignments)
        / len(assignments),
    }


def numeric(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    qc_rows = read_tsv(args.tree_qc)
    qc = {row["split_family"]: row for row in qc_rows}
    annotations: dict[str, dict[str, str]] = {}
    if args.annotations and args.annotations.is_file():
        annotations = {
            row["split_family"]: row for row in read_tsv(args.annotations)
            if row.get("split_family")
        }

    rows = []
    for family_dir in sorted(path for path in args.job_dir.iterdir() if path.is_dir()):
        family = family_dir.name
        family_qc = qc.get(family, {})
        annotation = annotations.get(family, {})
        null = best_run(family_dir, "null")
        alternative = best_run(family_dir, "alternative")
        complete = null is not None and alternative is not None
        warnings: list[str] = []

        if complete:
            null_omegas = null["omegas"]
            alt_omegas = alternative["omegas"]
            if len(null_omegas) != 2:
                warnings.append(f"null_expected_2_omegas_observed_{len(null_omegas)}")
            if len(alt_omegas) != 3:
                warnings.append(f"alternative_expected_3_omegas_observed_{len(alt_omegas)}")
            if len(null_omegas) == 2 and len(alt_omegas) == 3:
                omega_o = alt_omegas[0]
                omega_h = alt_omegas[1]
                omega_l = alt_omegas[2]
                omega_pooled_hl = null_omegas[1]
                ratio = omega_h / omega_l if omega_l > 0 else math.inf
            else:
                omega_o = omega_h = omega_l = omega_pooled_hl = ratio = None

            raw_delta = alternative["lnL"] - null["lnL"]
            if raw_delta < -0.01:
                warnings.append("alternative_lnL_below_null")
            lrt = max(0.0, 2.0 * raw_delta)
            p_value = chi1_sf(lrt)

            if null["n_successful_reps"] < 2:
                warnings.append("null_fewer_than_2_successful_reps")
            if alternative["n_successful_reps"] < 2:
                warnings.append("alternative_fewer_than_2_successful_reps")
            if null["n_reps_at_best_0.01"] < 2:
                warnings.append("null_best_optimum_not_reproduced")
            if alternative["n_reps_at_best_0.01"] < 2:
                warnings.append("alternative_best_optimum_not_reproduced")
            if null["lnL_rep_spread"] > 1:
                warnings.append("null_multistart_lnL_spread_gt_1")
            elif null["lnL_rep_spread"] > 0.01:
                warnings.append("null_multistart_lnL_spread_gt_0.01")
            if alternative["lnL_rep_spread"] > 1:
                warnings.append("alternative_multistart_lnL_spread_gt_1")
            elif alternative["lnL_rep_spread"] > 0.01:
                warnings.append("alternative_multistart_lnL_spread_gt_0.01")

            for label, omega in (
                ("omega_H", omega_h), ("omega_L", omega_l), ("omega_O", omega_o)
            ):
                if omega is not None and (omega <= 1e-4 or omega >= 100):
                    warnings.append(label + "_boundary")

            alignment_codons = numeric(family_qc.get("alignment_codons", ""))
            if alignment_codons is not None and alignment_codons < 100:
                warnings.append("short_alignment_lt_100_codons")
            if alternative["tree_dS"] is not None and alternative["tree_dS"] < 0.1:
                warnings.append("very_low_total_tree_dS")

            h_diag = class_diagnostics(alternative["branch_rows"], alt_omegas, 1)
            l_diag = class_diagnostics(alternative["branch_rows"], alt_omegas, 2)
            o_diag = class_diagnostics(alternative["branch_rows"], alt_omegas, 0)
            for label, diagnostic in (("HighT", h_diag), ("LowT", l_diag)):
                class_sum_ds = diagnostic["sum_dS"]
                if class_sum_ds is not None and class_sum_ds < 0.1:
                    warnings.append(label + "_class_sum_dS_lt_0.1")

            stable = (
                raw_delta >= -0.01
                and null["n_reps_at_best_0.01"] >= 2
                and alternative["n_reps_at_best_0.01"] >= 2
                and len(null_omegas) == 2
                and len(alt_omegas) == 3
            )
            convergence_status = "stable" if stable else "review"
        else:
            omega_o = omega_h = omega_l = omega_pooled_hl = ratio = None
            lrt = p_value = None
            h_diag = l_diag = o_diag = {
                "n_branches": 0, "sum_dN": None, "sum_dS": None,
                "median_dS": None, "zero_dS_fraction": None,
            }
            convergence_status = "incomplete"
            warnings.append("incomplete_models")

        rows.append({
            "split_family": family,
            "annotation": annotation.get(
                "Recommended_annotation", annotation.get("annotation", "")
            ),
            "KO": annotation.get("Integrated_KO", ""),
            "n_taxa": family_qc.get("n_taxa", ""),
            "alignment_codons": family_qc.get("alignment_codons", ""),
            "omega_HighT": "" if omega_h is None else omega_h,
            "omega_LowT": "" if omega_l is None else omega_l,
            "omega_Other": "" if omega_o is None else omega_o,
            "omega_pooled_HighT_LowT_null": (
                "" if omega_pooled_hl is None else omega_pooled_hl
            ),
            "omega_ratio_HighT_LowT": "" if ratio is None else ratio,
            "null_lnL": "" if null is None else null["lnL"],
            "alternative_lnL": "" if alternative is None else alternative["lnL"],
            "LRT_df": 1,
            "LRT": "" if lrt is None else lrt,
            "raw_P": "" if p_value is None else p_value,
            "BH_FDR": "",
            "direct_significant_FDR05": 0,
            "direction": (
                "" if ratio is None else
                ("HighT_lower" if ratio < 1 else "LowT_lower" if ratio > 1 else "equal")
            ),
            "convergence_status": convergence_status,
            "null_successful_reps": 0 if null is None else null["n_successful_reps"],
            "alternative_successful_reps": (
                0 if alternative is None else alternative["n_successful_reps"]
            ),
            "null_reps_at_best_0.01": (
                0 if null is None else null["n_reps_at_best_0.01"]
            ),
            "alternative_reps_at_best_0.01": (
                0 if alternative is None else alternative["n_reps_at_best_0.01"]
            ),
            "null_lnL_rep_spread": "" if null is None else null["lnL_rep_spread"],
            "alternative_lnL_rep_spread": (
                "" if alternative is None else alternative["lnL_rep_spread"]
            ),
            "null_all_rep_lnL": "" if null is None else null["all_rep_lnL"],
            "alternative_all_rep_lnL": (
                "" if alternative is None else alternative["all_rep_lnL"]
            ),
            "tree_length_dN": (
                "" if alternative is None or alternative["tree_dN"] is None
                else alternative["tree_dN"]
            ),
            "tree_length_dS": (
                "" if alternative is None or alternative["tree_dS"] is None
                else alternative["tree_dS"]
            ),
            "HighT_n_branches": h_diag["n_branches"],
            "HighT_sum_dN": "" if h_diag["sum_dN"] is None else h_diag["sum_dN"],
            "HighT_sum_dS": "" if h_diag["sum_dS"] is None else h_diag["sum_dS"],
            "HighT_median_branch_dS": (
                "" if h_diag["median_dS"] is None else h_diag["median_dS"]
            ),
            "HighT_zero_dS_fraction": (
                "" if h_diag["zero_dS_fraction"] is None
                else h_diag["zero_dS_fraction"]
            ),
            "LowT_n_branches": l_diag["n_branches"],
            "LowT_sum_dN": "" if l_diag["sum_dN"] is None else l_diag["sum_dN"],
            "LowT_sum_dS": "" if l_diag["sum_dS"] is None else l_diag["sum_dS"],
            "LowT_median_branch_dS": (
                "" if l_diag["median_dS"] is None else l_diag["median_dS"]
            ),
            "LowT_zero_dS_fraction": (
                "" if l_diag["zero_dS_fraction"] is None
                else l_diag["zero_dS_fraction"]
            ),
            "Other_n_branches": o_diag["n_branches"],
            "Other_sum_dN": "" if o_diag["sum_dN"] is None else o_diag["sum_dN"],
            "Other_sum_dS": "" if o_diag["sum_dS"] is None else o_diag["sum_dS"],
            "tree_qc_pass": family_qc.get("tree_qc_pass", ""),
            "warnings": ";".join(sorted(set(warnings))),
        })

    q_values = bh_adjust([
        float(row["raw_P"]) if row["raw_P"] != "" else None for row in rows
    ])
    for row, q_value in zip(rows, q_values):
        row["BH_FDR"] = "" if q_value is None else q_value
        row["direct_significant_FDR05"] = int(
            q_value is not None
            and q_value < 0.05
            and row["convergence_status"] == "stable"
        )

    fields = list(rows[0]) if rows else []
    write_tsv(args.outdir / "direct_HighT_LowT_omega_LRT.tsv", rows, fields)

    complete = sum(row["convergence_status"] != "incomplete" for row in rows)
    stable = sum(row["convergence_status"] == "stable" for row in rows)
    significant = sum(row["direct_significant_FDR05"] for row in rows)
    high_lower = sum(row["direction"] == "HighT_lower" for row in rows)
    low_lower = sum(row["direction"] == "LowT_lower" for row in rows)
    report = [
        "# Direct HighT-versus-LowT omega comparison",
        "",
        "## Model",
        "",
        "- Null: Other branches have omega_O; HighT and LowT within-crown branches share omega_HL.",
        "- Alternative: Other, HighT and LowT within-crown branches have omega_O, omega_H and omega_L.",
        "- The models differ by one free omega parameter; the LRT uses 1 df.",
        "- Crown scope includes terminal and internal descendant branches but excludes the entering stem.",
        "",
        "## Completion",
        "",
        f"- Families expected: {len(rows)}",
        f"- Families with both models complete: {complete}",
        f"- Families with reproducible multi-start optima: {stable}",
        f"- Direct per-gene tests significant at BH-FDR < 0.05: {significant}",
        f"- Direction among fitted genes: HighT lower = {high_lower}; LowT lower = {low_lower}",
        "",
        "## Interpretation boundary",
        "",
        "The existing separate-foreground ratio summarizes the global direction across genes. "
        "This nested direct model tests omega_H = omega_L for each gene while estimating a "
        "separate omega for all remaining branches. Only the direct LRT with BH correction "
        "supports gene-specific HighT-LowT significance.",
    ]
    (args.outdir / "DIRECT_OMEGA_SUMMARY.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print("\n".join(report[2:]))


if __name__ == "__main__":
    main()
