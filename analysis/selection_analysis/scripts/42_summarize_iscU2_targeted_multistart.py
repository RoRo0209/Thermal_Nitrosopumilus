#!/usr/bin/env python3
"""Summarize targeted iscU2 multi-start fits and conditionally update its LRT/FDR."""

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
KAPPA_RE = re.compile(r"kappa \(ts/tv\)\s*=\s*([0-9.Ee+-]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
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


def parse_omegas(text: str) -> list[float]:
    match = OMEGA_RE.search(text)
    if not match:
        return []
    values = []
    for token in match.group(1).split():
        try:
            values.append(float(token))
        except ValueError:
            pass
    return values


def chi1_sf(statistic: float) -> float:
    return math.erfc(math.sqrt(max(statistic, 0.0) / 2.0))


def bh_adjust(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    adjusted = [1.0] * len(values)
    running = 1.0
    for reverse_index in range(len(order) - 1, -1, -1):
        index = order[reverse_index]
        rank = reverse_index + 1
        running = min(running, values[index] * len(values) / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    manifest = read_tsv(args.root / "iscU2_multistart_manifest.tsv")
    original_path = (
        args.root / "inputs" / "direct_HighT_LowT_omega_LRT_original.tsv"
    )
    original = read_tsv(original_path)
    old = next(row for row in original if row["split_family"] == "iscU2")

    runs = []
    for item in manifest:
        run_dir = args.root / "jobs" / "iscU2" / item["run_id"]
        result = run_dir / "results.txt"
        row = dict(item)
        row.update(
            {
                "completed": 0,
                "final_lnL": "",
                "final_kappa": "",
                "final_omega_Other": "",
                "final_omega_HighT_or_pooled": "",
                "final_omega_LowT": "",
                "final_HighT_LowT_ratio": "",
                "direction": "",
                "within_0.01_of_model_best": 0,
            }
        )
        if result.is_file() and result.stat().st_size:
            text = result.read_text(encoding="utf-8", errors="replace")
            lnl = LNL_RE.search(text)
            kappa = KAPPA_RE.search(text)
            omegas = parse_omegas(text)
            if "Time used:" in text and lnl:
                row["completed"] = 1
                row["final_lnL"] = float(lnl.group(1))
                row["final_kappa"] = float(kappa.group(1)) if kappa else ""
                if item["model"] == "null" and len(omegas) == 2:
                    row["final_omega_Other"] = omegas[0]
                    row["final_omega_HighT_or_pooled"] = omegas[1]
                if item["model"] == "alternative" and len(omegas) == 3:
                    row["final_omega_Other"] = omegas[0]
                    row["final_omega_HighT_or_pooled"] = omegas[1]
                    row["final_omega_LowT"] = omegas[2]
                    if omegas[2] > 0:
                        ratio = omegas[1] / omegas[2]
                        row["final_HighT_LowT_ratio"] = ratio
                        row["direction"] = (
                            "HighT_lower" if ratio < 1
                            else "LowT_lower" if ratio > 1
                            else "equal"
                        )
        runs.append(row)

    completed = [row for row in runs if row["completed"] == 1]
    if len(completed) != len(runs):
        raise RuntimeError(
            f"Expected {len(runs)} completed runs, observed {len(completed)}"
        )

    summaries = {}
    for model in ("null", "alternative"):
        model_runs = [row for row in completed if row["model"] == model]
        best = max(float(row["final_lnL"]) for row in model_runs)
        at_best = []
        for row in model_runs:
            flag = best - float(row["final_lnL"]) <= 0.01
            row["within_0.01_of_model_best"] = int(flag)
            if flag:
                at_best.append(row)
        summaries[model] = {"best": best, "at_best": at_best}

    null_best = summaries["null"]["best"]
    alt_best = summaries["alternative"]["best"]
    old_null = float(old["null_lnL"])
    old_alt = float(old["alternative_lnL"])
    best_changed = (
        abs(null_best - old_null) > 0.01 or abs(alt_best - old_alt) > 0.01
    )

    alt_cluster = summaries["alternative"]["at_best"]
    ratios = [
        float(row["final_HighT_LowT_ratio"])
        for row in alt_cluster
        if row["final_HighT_LowT_ratio"] != ""
    ]
    directions = {
        row["direction"] for row in alt_cluster if row["direction"]
    }
    omega_h = [
        float(row["final_omega_HighT_or_pooled"])
        for row in alt_cluster
        if row["final_omega_HighT_or_pooled"] != ""
    ]
    omega_l = [
        float(row["final_omega_LowT"])
        for row in alt_cluster
        if row["final_omega_LowT"] != ""
    ]

    if best_changed:
        lrt = max(0.0, 2.0 * (alt_best - null_best))
        p_value = chi1_sf(lrt)
        updated = [dict(row) for row in original]
        for row in updated:
            if row["split_family"] == "iscU2":
                row["null_lnL"] = null_best
                row["alternative_lnL"] = alt_best
                row["LRT"] = lrt
                row["raw_P"] = p_value
        q_values = bh_adjust([float(row["raw_P"]) for row in updated])
        for row, q_value in zip(updated, q_values):
            row["BH_FDR"] = q_value
        updated_isc = next(
            row for row in updated if row["split_family"] == "iscU2"
        )
        fdr = float(updated_isc["BH_FDR"])
        write_tsv(
            args.outdir / "direct_HighT_LowT_omega_LRT_updated_if_needed.tsv",
            updated,
            list(updated[0]),
        )
        fdr_action = "recomputed because an accepted best lnL changed"
    else:
        lrt = float(old["LRT"])
        p_value = float(old["raw_P"])
        fdr = float(old["BH_FDR"])
        fdr_action = "retained because both accepted best lnL values were unchanged"

    null_replicates = len(summaries["null"]["at_best"])
    alt_replicates = len(summaries["alternative"]["at_best"])
    stable = (
        null_replicates >= 2
        and alt_replicates >= 2
        and len(directions) == 1
        and bool(ratios)
    )
    preferred = null_replicates >= 3 and alt_replicates >= 3

    write_tsv(
        args.outdir / "iscU2_all_targeted_multistart_runs.tsv",
        runs,
        list(runs[0]),
    )
    report = f"""# iscU2 targeted multi-start convergence

## Grid

- Null runs: 30
- Alternative runs: 30
- Initial omega values: 0.005, 0.02, 0.05, 0.2, 1, 5
- Initial kappa values: 0.25, 0.75, 2, 5, 12

## Reproducibility

- Null best lnL: {null_best:.6f}; runs within 0.01: {null_replicates}/30
- Alternative best lnL: {alt_best:.6f}; runs within 0.01: {alt_replicates}/30
- Independently reproduced at least twice: {"yes" if stable else "no"}
- Independently reproduced at least three times for both models: {"yes" if preferred else "no"}
- Previous null best lnL: {old_null:.6f}
- Previous alternative best lnL: {old_alt:.6f}
- Accepted best likelihood changed by >0.01: {"yes" if best_changed else "no"}

## Alternative-model effect stability

- Direction(s) among best-cluster runs: {", ".join(sorted(directions)) or "missing"}
- omega HighT range among best-cluster runs: {min(omega_h):.6g} to {max(omega_h):.6g}
- omega LowT range among best-cluster runs: {min(omega_l):.6g} to {max(omega_l):.6g}
- HighT/LowT ratio range: {min(ratios):.6g} to {max(ratios):.6g}

## Test

- LRT (1 df): {lrt:.6f}
- Raw P: {p_value:.8g}
- BH-FDR: {fdr:.8g}
- FDR action: {fdr_action}
- Final convergence decision: {"accepted" if stable else "remains unstable"}
"""
    (args.outdir / "ISCU2_TARGETED_CONVERGENCE.md").write_text(
        report, encoding="utf-8", newline="\n"
    )
    print(report)


if __name__ == "__main__":
    main()
