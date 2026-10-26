#!/usr/bin/env python3
"""Summarize replicated codeml branch/branch-site runs, LRTs, BEB sites and annotations."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


LNL_RE = re.compile(r"lnL\(.*?\):\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)")
OMEGA_RE = re.compile(r"w \(dN/dS\) for branches:\s*(.+)")
BEB_RE = re.compile(r"^\s*(\d+)\s+([A-Z*?-])\s+([01](?:\.\d+)?)\s*(\*{0,2})", re.MULTILINE)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--job-dir", type=Path, required=True)
    p.add_argument("--family-annotations", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    return p.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def best_run(ecotype_dir: Path, prefix: str):
    runs = []
    for run_dir in sorted(ecotype_dir.glob(prefix + "_rep*")):
        result = run_dir / "results.txt"
        if not result.is_file():
            continue
        text = result.read_text(encoding="utf-8", errors="replace")
        match = LNL_RE.search(text)
        if not match:
            continue
        runs.append((float(match.group(1)), run_dir, text))
    if not runs:
        return None
    runs.sort(key=lambda item: item[0], reverse=True)
    best = runs[0]
    return {
        "lnL": best[0],
        "run_dir": str(best[1]),
        "text": best[2],
        "n_successful_reps": len(runs),
        "lnL_rep_spread": max(x[0] for x in runs) - min(x[0] for x in runs),
    }


def chi1_sf(statistic: float) -> float:
    return math.erfc(math.sqrt(max(statistic, 0.0) / 2.0))


def bh_adjust(values: list[float | None]) -> list[float | None]:
    indexed = [(i, p) for i, p in enumerate(values) if p is not None and math.isfinite(p)]
    indexed.sort(key=lambda item: item[1])
    adjusted: dict[int, float] = {}
    running = 1.0
    m = len(indexed)
    for rank_from_end in range(m - 1, -1, -1):
        index, p = indexed[rank_from_end]
        rank = rank_from_end + 1
        running = min(running, p * m / rank)
        adjusted[index] = min(running, 1.0)
    return [adjusted.get(i) for i in range(len(values))]


def branch_omegas(text: str):
    match = OMEGA_RE.search(text)
    if not match:
        return None, None
    values = []
    for token in match.group(1).split():
        try:
            values.append(float(token))
        except ValueError:
            pass
    if len(values) < 2:
        return None, None
    return values[0], values[-1]


def beb_sites(text: str):
    marker = text.rfind("Bayes Empirical Bayes (BEB) analysis")
    if marker < 0:
        return []
    section = text[marker:]
    return [
        {"codon_position": int(m.group(1)), "amino_acid": m.group(2),
         "posterior_probability": float(m.group(3)), "stars": m.group(4)}
        for m in BEB_RE.finditer(section)
    ]


def write_tsv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    annotation_rows = read_tsv(args.family_annotations)
    annotations = {row["split_family"]: row for row in annotation_rows}
    file_to_family = {}
    for row in annotation_rows:
        family = row["split_family"]
        file_to_family[row.get("file_stem") or re.sub(r"[^A-Za-z0-9_.-]+", "_", family)] = family

    rows = []
    beb_rows = []
    for family_dir in sorted(path for path in args.job_dir.iterdir() if path.is_dir()):
        stem = family_dir.name
        family = file_to_family.get(stem, stem)
        annotation = annotations.get(family, {})
        for ecotype in ("HighT", "LowT"):
            ecotype_dir = family_dir / ecotype
            branch_null = best_run(ecotype_dir, "branch_null")
            branch_alt = best_run(ecotype_dir, "branch_alt")
            bs_null = best_run(ecotype_dir, "branchsite_null")
            bs_alt = best_run(ecotype_dir, "branchsite_alt")
            complete = all(x is not None for x in (branch_null, branch_alt, bs_null, bs_alt))
            if complete:
                branch_lrt = max(0.0, 2 * (branch_alt["lnL"] - branch_null["lnL"]))
                branch_p = chi1_sf(branch_lrt)
                bs_lrt = max(0.0, 2 * (bs_alt["lnL"] - bs_null["lnL"]))
                bs_p_chi1 = chi1_sf(bs_lrt)
                bs_p_mixture = 0.5 * bs_p_chi1
                omega_background, omega_foreground = branch_omegas(branch_alt["text"])
                sites = beb_sites(bs_alt["text"])
                for site in sites:
                    if site["posterior_probability"] >= 0.50:
                        beb_rows.append({
                            "split_family": family,
                            "file_stem": stem,
                            "ecotype": ecotype,
                            **site,
                            "annotation": annotation.get("annotation", ""),
                            "non_unique_gene_name": annotation.get("non_unique_gene_name", ""),
                        })
            else:
                branch_lrt = branch_p = bs_lrt = bs_p_chi1 = bs_p_mixture = None
                omega_background = omega_foreground = None
                sites = []
            rows.append({
                "split_family": family,
                "file_stem": stem,
                "ecotype": ecotype,
                "all_models_complete": int(complete),
                "branch_null_lnL": branch_null["lnL"] if branch_null else "",
                "branch_alt_lnL": branch_alt["lnL"] if branch_alt else "",
                "branch_LRT": branch_lrt if branch_lrt is not None else "",
                "branch_p_chi1": branch_p if branch_p is not None else "",
                "omega_background": omega_background if omega_background is not None else "",
                "omega_foreground": omega_foreground if omega_foreground is not None else "",
                "branch_alt_successful_reps": branch_alt["n_successful_reps"] if branch_alt else 0,
                "branch_alt_lnL_rep_spread": branch_alt["lnL_rep_spread"] if branch_alt else "",
                "branchsite_null_lnL": bs_null["lnL"] if bs_null else "",
                "branchsite_alt_lnL": bs_alt["lnL"] if bs_alt else "",
                "branchsite_LRT": bs_lrt if bs_lrt is not None else "",
                "branchsite_p_chi1": bs_p_chi1 if bs_p_chi1 is not None else "",
                "branchsite_p_mixture": bs_p_mixture if bs_p_mixture is not None else "",
                "branchsite_alt_successful_reps": bs_alt["n_successful_reps"] if bs_alt else 0,
                "branchsite_alt_lnL_rep_spread": bs_alt["lnL_rep_spread"] if bs_alt else "",
                "n_BEB_ge_0.95": sum(site["posterior_probability"] >= 0.95 for site in sites),
                "n_BEB_ge_0.99": sum(site["posterior_probability"] >= 0.99 for site in sites),
                "non_unique_gene_name": annotation.get("non_unique_gene_name", ""),
                "annotation": annotation.get("annotation", ""),
                "nosplit_ALE_families": annotation.get("nosplit_ALE_families", ""),
            })

    for ecotype in ("HighT", "LowT"):
        indexes = [i for i, row in enumerate(rows) if row["ecotype"] == ecotype]
        branch_q = bh_adjust([
            float(rows[i]["branch_p_chi1"]) if rows[i]["branch_p_chi1"] != "" else None
            for i in indexes
        ])
        bs_q = bh_adjust([
            float(rows[i]["branchsite_p_mixture"]) if rows[i]["branchsite_p_mixture"] != "" else None
            for i in indexes
        ])
        for local, global_i in enumerate(indexes):
            rows[global_i]["branch_FDR"] = branch_q[local] if branch_q[local] is not None else ""
            rows[global_i]["branchsite_FDR_mixture"] = bs_q[local] if bs_q[local] is not None else ""
            rows[global_i]["branch_significant_FDR05"] = int(branch_q[local] is not None and branch_q[local] < 0.05)
            rows[global_i]["branchsite_significant_FDR05"] = int(bs_q[local] is not None and bs_q[local] < 0.05)

    fields = [
        "split_family", "file_stem", "ecotype", "all_models_complete",
        "branch_null_lnL", "branch_alt_lnL", "branch_LRT", "branch_p_chi1", "branch_FDR",
        "omega_background", "omega_foreground", "branch_alt_successful_reps", "branch_alt_lnL_rep_spread",
        "branch_significant_FDR05", "branchsite_null_lnL", "branchsite_alt_lnL", "branchsite_LRT",
        "branchsite_p_chi1", "branchsite_p_mixture", "branchsite_FDR_mixture",
        "branchsite_alt_successful_reps", "branchsite_alt_lnL_rep_spread",
        "branchsite_significant_FDR05", "n_BEB_ge_0.95", "n_BEB_ge_0.99",
        "non_unique_gene_name", "annotation", "nosplit_ALE_families",
    ]
    write_tsv(args.outdir / "codeml_summary_all.tsv", rows, fields)
    beb_fields = [
        "split_family", "file_stem", "ecotype", "codon_position", "amino_acid",
        "posterior_probability", "stars", "non_unique_gene_name", "annotation",
    ]
    write_tsv(args.outdir / "branchsite_BEB_sites.tsv", beb_rows, beb_fields)

    complete = sum(int(row["all_models_complete"]) for row in rows)
    branch_sig = sum(int(row.get("branch_significant_FDR05", 0)) for row in rows)
    bs_sig = sum(int(row.get("branchsite_significant_FDR05", 0)) for row in rows)
    report = [
        "# codeml branch and branch-site summary",
        "",
        f"- Family-ecotype tests expected: {len(rows)}",
        f"- Family-ecotype tests with all models complete: {complete}",
        f"- Branch-model tests significant at BH-FDR < 0.05: {branch_sig}",
        f"- Branch-site tests significant at mixture-test BH-FDR < 0.05: {bs_sig}",
        f"- BEB sites reported at posterior probability >= 0.50: {len(beb_rows)}",
        "",
        "Branch tests use a one-ratio versus two-ratio LRT (1 df). Branch-site tests report both the "
        "conservative chi-square(1) p value and the standard 50:50 point-mass/chi-square(1) mixture p value. "
        "FDR is controlled separately for HighT and LowT.",
    ]
    (args.outdir / "CODEML_RESULTS_SUMMARY.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report[2:]))


if __name__ == "__main__":
    main()
