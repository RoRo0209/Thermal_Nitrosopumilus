#!/usr/bin/env python3
"""Audit and summarize the downloaded HighT/LowT codeml results."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def f(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in ("", None) else math.nan


def i(row: dict[str, str], key: str) -> int:
    value = row.get(key, "")
    return int(float(value)) if value not in ("", None) else 0


def bh_adjust(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=lambda idx: values[idx])
    adjusted = [1.0] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        idx = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[idx] * n / rank)
        adjusted[idx] = min(1.0, running)
    return adjusted


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float) -> str:
    if math.isnan(value):
        return "NA"
    if value == 0:
        return "0"
    if abs(value) < 0.001:
        return f"{value:.2e}"
    return f"{value:.4g}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary_path = args.input_dir / "selection_analysis/results/09_codeml_summary_soft85/codeml_summary_all.tsv"
    beb_path = args.input_dir / "selection_analysis/results/09_codeml_summary_soft85/branchsite_BEB_sites.tsv"
    qc_path = args.input_dir / "selection_analysis/results/08_codeml_jobs_soft85/codon_alignment_qc.tsv"
    summary = read_tsv(summary_path)
    beb = read_tsv(beb_path)
    qc = read_tsv(qc_path)

    keys = [(r["split_family"], r["ecotype"]) for r in summary]
    key_counts = Counter(keys)
    duplicates = [key for key, count in key_counts.items() if count > 1]
    families = sorted({r["split_family"] for r in summary})
    family_ecotypes = defaultdict(set)
    for row in summary:
        family_ecotypes[row["split_family"]].add(row["ecotype"])
    incomplete_pairs = [fam for fam, ecotypes in family_ecotypes.items() if ecotypes != {"HighT", "LowT"}]

    fdr_mismatches: list[str] = []
    for ecotype in ("HighT", "LowT"):
        subset = [r for r in summary if r["ecotype"] == ecotype]
        for p_col, q_col in (("branch_p_chi1", "branch_FDR"), ("branchsite_p_mixture", "branchsite_FDR_mixture")):
            recalculated = bh_adjust([f(r, p_col) for r in subset])
            for row, q_new in zip(subset, recalculated):
                if abs(q_new - f(row, q_col)) > 1e-9:
                    fdr_mismatches.append(f"{row['split_family']}:{ecotype}:{q_col}")

    qc_eligible = [r for r in qc if i(r, "codeml_eligible") == 1]
    qc_ineligible = [r for r in qc if i(r, "codeml_eligible") == 0]
    all_complete = sum(i(r, "all_models_complete") == 1 for r in summary)
    alt_rep_fail = [r for r in summary if i(r, "branch_alt_successful_reps") < 3 or i(r, "branchsite_alt_successful_reps") < 3]
    branch_spread = [r for r in summary if f(r, "branch_alt_lnL_rep_spread") > 0.01]
    branchsite_spread = [r for r in summary if f(r, "branchsite_alt_lnL_rep_spread") > 0.01]

    branch_sig = [r for r in summary if i(r, "branch_significant_FDR05") == 1]
    branchsite_sig = [r for r in summary if i(r, "branchsite_significant_FDR05") == 1]

    def family_sets(rows: list[dict[str, str]]) -> dict[str, set[str]]:
        return {
            eco: {r["split_family"] for r in rows if r["ecotype"] == eco}
            for eco in ("HighT", "LowT")
        }

    branch_sets = family_sets(branch_sig)
    branchsite_sets = family_sets(branchsite_sig)

    beb_by_key = defaultdict(list)
    for row in beb:
        beb_by_key[(row["split_family"], row["ecotype"])].append(row)

    candidate_rows: list[dict[str, object]] = []
    for row in sorted(branchsite_sig, key=lambda r: (f(r, "branchsite_FDR_mixture"), r["ecotype"], r["split_family"])):
        sites = beb_by_key[(row["split_family"], row["ecotype"])]
        sites95 = [s for s in sites if f(s, "posterior_probability") >= 0.95]
        sites99 = [s for s in sites if f(s, "posterior_probability") >= 0.99]
        candidate_rows.append({
            "split_family": row["split_family"],
            "ecotype": row["ecotype"],
            "branchsite_LRT": row["branchsite_LRT"],
            "branchsite_p_mixture": row["branchsite_p_mixture"],
            "branchsite_FDR_mixture": row["branchsite_FDR_mixture"],
            "n_BEB_ge_0.95": len(sites95),
            "n_BEB_ge_0.99": len(sites99),
            "BEB_sites_ge_0.95": ";".join(f"{s['amino_acid']}{s['codon_position']}({float(s['posterior_probability']):.3f})" for s in sites95),
            "annotation": row["annotation"],
            "nosplit_ALE_families": row["nosplit_ALE_families"],
        })

    branch_rows: list[dict[str, object]] = []
    for row in sorted(branch_sig, key=lambda r: (f(r, "branch_FDR"), r["ecotype"], r["split_family"])):
        bg = f(row, "omega_background")
        fg = f(row, "omega_foreground")
        direction = "higher" if fg > bg else "lower" if fg < bg else "equal"
        branch_rows.append({
            "split_family": row["split_family"],
            "ecotype": row["ecotype"],
            "branch_FDR": row["branch_FDR"],
            "omega_background": row["omega_background"],
            "omega_foreground": row["omega_foreground"],
            "foreground_vs_background": direction,
            "annotation": row["annotation"],
            "nosplit_ALE_families": row["nosplit_ALE_families"],
        })

    write_tsv(
        args.output_dir / "branchsite_significant_candidates.tsv",
        candidate_rows,
        ["split_family", "ecotype", "branchsite_LRT", "branchsite_p_mixture", "branchsite_FDR_mixture", "n_BEB_ge_0.95", "n_BEB_ge_0.99", "BEB_sites_ge_0.95", "annotation", "nosplit_ALE_families"],
    )
    write_tsv(
        args.output_dir / "branch_significant_rate_shifts.tsv",
        branch_rows,
        ["split_family", "ecotype", "branch_FDR", "omega_background", "omega_foreground", "foreground_vs_background", "annotation", "nosplit_ALE_families"],
    )

    high_conf = [r for r in candidate_rows if int(r["n_BEB_ge_0.95"]) > 0]
    strong_conf = [r for r in candidate_rows if int(r["n_BEB_ge_0.99"]) > 0]
    branch_direction = Counter((r["ecotype"], r["foreground_vs_background"]) for r in branch_rows)
    unannotated_site = [r for r in candidate_rows if not str(r["annotation"]).strip() or "hypothetical" in str(r["annotation"]).lower()]

    top_lines = []
    for row in candidate_rows[:15]:
        top_lines.append(
            f"| {row['split_family']} | {row['ecotype']} | {fmt(float(row['branchsite_FDR_mixture']))} | "
            f"{row['n_BEB_ge_0.95']} | {row['n_BEB_ge_0.99']} | {row['BEB_sites_ge_0.95'] or '—'} | {row['annotation'] or 'unannotated'} |"
        )

    report = f"""# HighT/LowT codeml结果核查与初步统计

## 数据与粒度

- codeml汇总：{len(summary)}行，代表{len(families)}个家族 × 2个生态型前景检验。
- BEB记录：{len(beb)}个位点（原始表保留阈值为posterior probability ≥ 0.50）。
- 密码子比对QC：{len(qc)}个家族，其中{len(qc_eligible)}个进入codeml，{len(qc_ineligible)}个未进入。

## 完整性与质量核查

- all_models_complete：{all_complete}/{len(summary)}。
- 重复的family–ecotype键：{len(duplicates)}。
- 缺少HighT或LowT配对的家族：{len(incomplete_pairs)}。
- 少于3个alternative起始值成功的检验：{len(alt_rep_fail)}。
- 重新计算BH-FDR后不一致的记录：{len(fdr_mismatches)}。
- branch alternative重复间lnL极差 > 0.01：{len(branch_spread)}/{len(summary)}。
- branch-site alternative重复间lnL极差 > 0.01：{len(branchsite_spread)}/{len(summary)}。

## Branch model：前景与背景的整体ω不同

- FDR < 0.05：{len(branch_sig)}/{len(summary)}个family–ecotype检验，涉及{len(branch_sets['HighT'] | branch_sets['LowT'])}个不同家族。
- HighT：{len(branch_sets['HighT'])}个；其中前景ω降低{branch_direction[('HighT','lower')]}个、升高{branch_direction[('HighT','higher')]}个。
- LowT：{len(branch_sets['LowT'])}个；其中前景ω降低{branch_direction[('LowT','lower')]}个、升高{branch_direction[('LowT','higher')]}个。
- 两个生态型均显著：{len(branch_sets['HighT'] & branch_sets['LowT'])}个家族；仅HighT：{len(branch_sets['HighT'] - branch_sets['LowT'])}个；仅LowT：{len(branch_sets['LowT'] - branch_sets['HighT'])}个。

这里的显著性只说明指定前景枝与其余所有背景枝的ω不同：HighT检验的背景包含LowT和Other NP，LowT检验的背景包含HighT和Other NP。前景ω降低通常与更强的净化选择一致；前景ω升高可由选择约束放松或正选择造成，不能仅凭branch model区分。

## Branch-site model：前景stem上少数密码子正选择

- mixture-test FDR < 0.05：{len(branchsite_sig)}/{len(summary)}个family–ecotype检验，涉及{len(branchsite_sets['HighT'] | branchsite_sets['LowT'])}个不同家族。
- HighT：{len(branchsite_sets['HighT'])}个；LowT：{len(branchsite_sets['LowT'])}个。
- 两个生态型均显著：{len(branchsite_sets['HighT'] & branchsite_sets['LowT'])}个家族；仅HighT：{len(branchsite_sets['HighT'] - branchsite_sets['LowT'])}个；仅LowT：{len(branchsite_sets['LowT'] - branchsite_sets['HighT'])}个。
- 同时至少有1个BEB ≥ 0.95位点：{len(high_conf)}/{len(candidate_rows)}个显著检验。
- 同时至少有1个BEB ≥ 0.99位点：{len(strong_conf)}/{len(candidate_rows)}个显著检验。
- 显著候选中未注释/假定蛋白：{len(unannotated_site)}/{len(candidate_rows)}。

## 按FDR排序的前15个branch-site检验

| Family | Ecotype | FDR | BEB≥0.95 | BEB≥0.99 | 高可信位点 | Annotation |
|---|---|---:|---:|---:|---|---|
{chr(10).join(top_lines)}

## 解释边界与下一步

1. HighT和LowT分别作为前景，其余所有枝作为背景；两套独立检验不能直接证明HighT与LowT彼此不同。
2. 正文候选应优先要求branch-site FDR < 0.05且至少一个BEB ≥ 0.95；BEB ≥ 0.50只用于探索，不宜作为强证据。
3. lnL重复间极差 > 0.01的模型需要检查是否存在不同局部最优值；最终候选尤其应优先核查。
4. 下一步应将候选与群体pN/pS、ALE事件和功能注释整合，再决定主图和候选基因树。
"""
    (args.output_dir / "CODEML_ANALYSIS_REPORT_CN.md").write_text(report, encoding="utf-8")

    print(report)
    if qc_ineligible:
        print("\nIneligible QC families:", ", ".join(r["split_family"] for r in qc_ineligible))
    if branch_spread:
        print("Branch lnL spread >0.01 (top):")
        for row in sorted(branch_spread, key=lambda r: f(r, "branch_alt_lnL_rep_spread"), reverse=True)[:10]:
            print(row["split_family"], row["ecotype"], row["branch_alt_lnL_rep_spread"])
    if branchsite_spread:
        print("Branch-site lnL spread >0.01 (top):")
        for row in sorted(branchsite_spread, key=lambda r: f(r, "branchsite_alt_lnL_rep_spread"), reverse=True)[:10]:
            print(row["split_family"], row["ecotype"], row["branchsite_alt_lnL_rep_spread"])


if __name__ == "__main__":
    main()
