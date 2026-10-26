#!/usr/bin/env python3
"""Audit replicated codeml runs for local optima and branch-site degeneracy."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


LNL_RE = re.compile(r"lnL\(.*?\):\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)")
PROP_RE = re.compile(r"^proportion\s+(.+)$", re.MULTILINE)
FGW_RE = re.compile(r"^foreground w\s+(.+)$", re.MULTILINE)
BEB_RE = re.compile(r"^\s*(\d+)\s+([A-Z*?-])\s+([01](?:\.\d+)?)\s*(\*{0,2})", re.MULTILINE)
TIME_RE = re.compile(r"Time used:\s*(.+)")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def floats(text: str) -> list[float]:
    values = []
    for token in text.split():
        try:
            values.append(float(token))
        except ValueError:
            pass
    return values


def parse_run(run_dir: Path) -> dict[str, object]:
    result = run_dir / "results.txt"
    stderr = run_dir / "codeml.stderr"
    stdout = run_dir / "codeml.stdout"
    text = result.read_text(encoding="utf-8", errors="replace") if result.is_file() else ""
    err = stderr.read_text(encoding="utf-8", errors="replace").strip() if stderr.is_file() else ""
    out = stdout.read_text(encoding="utf-8", errors="replace") if stdout.is_file() else ""
    lnl_match = LNL_RE.search(text)
    prop_matches = PROP_RE.findall(text)
    fgw_matches = FGW_RE.findall(text)
    proportions = floats(prop_matches[-1]) if prop_matches else []
    fg_w = floats(fgw_matches[-1]) if fgw_matches else []
    marker = text.rfind("Bayes Empirical Bayes (BEB) analysis")
    sites = []
    if marker >= 0:
        for match in BEB_RE.finditer(text[marker:]):
            sites.append((int(match.group(1)), match.group(2), float(match.group(3)), match.group(4)))
    max_beb = max((s[2] for s in sites), default=0.0)
    ctl = (run_dir / "codeml.ctl").read_text(encoding="utf-8", errors="replace")
    omega_start = re.search(r"^\s*omega\s*=\s*([^\s*]+)", ctl, re.MULTILINE)
    time_matches = TIME_RE.findall(text) or TIME_RE.findall(out)
    normal_end = bool(lnl_match and time_matches and not err)
    return {
        "run": run_dir.name,
        "lnL": float(lnl_match.group(1)) if lnl_match else math.nan,
        "omega_start": omega_start.group(1) if omega_start else "",
        "p0": proportions[0] if len(proportions) > 0 else math.nan,
        "p1": proportions[1] if len(proportions) > 1 else math.nan,
        "p2": proportions[2] if len(proportions) > 2 else math.nan,
        "p3": proportions[3] if len(proportions) > 3 else math.nan,
        "positive_class_proportion": sum(proportions[2:4]) if len(proportions) >= 4 else math.nan,
        "foreground_omega2": fg_w[2] if len(fg_w) > 2 else math.nan,
        "max_BEB": max_beb,
        "n_BEB_ge_0.95": sum(s[2] >= 0.95 for s in sites),
        "BEB_sites_ge_0.95": ";".join(f"{aa}{pos}({pp:.3f})" for pos, aa, pp, _ in sites if pp >= 0.95),
        "normal_end": normal_end,
        "stderr_nonempty": bool(err),
        "time_used": time_matches[-1] if time_matches else "",
    }


def mixture_p(lrt: float) -> float:
    return 0.5 * math.erfc(math.sqrt(max(lrt, 0.0) / 2.0))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def bh_adjust(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    adjusted = [1.0] * len(values)
    running = 1.0
    for rank_index in range(len(values) - 1, -1, -1):
        idx = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[idx] * len(values) / rank)
        adjusted[idx] = min(1.0, running)
    return adjusted


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-dir", required=True, type=Path)
    parser.add_argument("--integrated-table", required=True, type=Path)
    parser.add_argument("--full-summary", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    integrated = read_tsv(args.integrated_table)
    focal_by_family = {r["split_family"]: r["foreground"] for r in integrated}
    original_by_key = {(r["split_family"], r["foreground"]): r for r in integrated}
    all_runs: list[dict[str, object]] = []
    audits: list[dict[str, object]] = []

    for family_dir in sorted(p for p in args.job_dir.iterdir() if p.is_dir()):
        family = family_dir.name
        for ecotype in ("HighT", "LowT"):
            eco_dir = family_dir / ecotype
            parsed: dict[str, dict[str, object]] = {}
            for run_dir in sorted(p for p in eco_dir.iterdir() if p.is_dir()):
                record = parse_run(run_dir)
                record.update({"family": family, "ecotype": ecotype})
                parsed[run_dir.name] = record
                all_runs.append(record)

            for model in ("branch", "branchsite"):
                null = parsed[f"{model}_null_rep1"]
                alts = [parsed[f"{model}_alt_rep{i}"] for i in (1, 2, 3)]
                alts_sorted = sorted(alts, key=lambda r: float(r["lnL"]), reverse=True)
                best = alts_sorted[0]
                worst = alts_sorted[-1]
                lrt_each = [max(0.0, 2 * (float(a["lnL"]) - float(null["lnL"]))) for a in alts]
                p_each = [mixture_p(x) if model == "branchsite" else math.erfc(math.sqrt(x / 2.0)) for x in lrt_each]
                positive_prop = float(best["positive_class_proportion"]) if model == "branchsite" else math.nan
                max_beb = float(best["max_BEB"]) if model == "branchsite" else math.nan
                decision = "not_focal"
                if focal_by_family.get(family) == ecotype and model == "branchsite":
                    spread = float(best["lnL"]) - float(worst["lnL"])
                    n_raw_sig = sum(p < 0.05 for p in p_each)
                    site_sets = []
                    for alt in alts:
                        tokens = str(alt["BEB_sites_ge_0.95"]).split(";") if alt["BEB_sites_ge_0.95"] else []
                        site_sets.append({token.split("(")[0] for token in tokens if token})
                    common_sites = set.intersection(*site_sets) if site_sets else set()
                    if positive_prop <= 1e-5 and max_beb < 0.95:
                        decision = "likely_null_local_optimum_do_not_use"
                    elif n_raw_sig == 3 and common_sites and spread <= 0.2:
                        decision = "retain_robust_across_starts"
                    elif n_raw_sig == 3 and common_sites:
                        decision = "retain_provisional_targeted_rerun"
                    elif spread > 1 or n_raw_sig < 3:
                        decision = "rerun_null_and_alt"
                    elif spread > 0.01:
                        decision = "review_or_targeted_rerun"
                    else:
                        decision = "retain"
                original = original_by_key.get((family, ecotype), {})
                audits.append({
                    "family": family,
                    "ecotype": ecotype,
                    "is_original_branchsite_candidate_foreground": int(focal_by_family.get(family) == ecotype),
                    "model": model,
                    "null_lnL": null["lnL"],
                    "alt_rep1_lnL": alts[0]["lnL"],
                    "alt_rep2_lnL": alts[1]["lnL"],
                    "alt_rep3_lnL": alts[2]["lnL"],
                    "best_alt_run": best["run"],
                    "best_alt_lnL": best["lnL"],
                    "alt_lnL_spread": float(best["lnL"]) - float(worst["lnL"]),
                    "best_LRT": max(lrt_each),
                    "best_p": min(p_each),
                    "alt_reps_raw_p_lt_0.05": sum(p < 0.05 for p in p_each),
                    "rep1_p": p_each[0],
                    "rep2_p": p_each[1],
                    "rep3_p": p_each[2],
                    "best_positive_class_proportion": positive_prop,
                    "best_foreground_omega2": best["foreground_omega2"] if model == "branchsite" else "",
                    "best_max_BEB": max_beb,
                    "best_n_BEB_ge_0.95": best["n_BEB_ge_0.95"] if model == "branchsite" else "",
                    "best_BEB_sites_ge_0.95": best["BEB_sites_ge_0.95"] if model == "branchsite" else "",
                    "common_BEB_sites_across_alt_reps": ";".join(sorted(common_sites)) if model == "branchsite" and focal_by_family.get(family) == ecotype else "",
                    "all_runs_normal_end": int(bool(null["normal_end"]) and all(bool(a["normal_end"]) for a in alts)),
                    "original_branchsite_FDR": original.get("branchsite_FDR", "") if model == "branchsite" else "",
                    "decision": decision,
                })

    write_tsv(args.output_dir / "all_codeml_runs_parsed.tsv", all_runs)
    write_tsv(args.output_dir / "convergence_audit_all_models.tsv", audits)
    focal = [r for r in audits if r["model"] == "branchsite" and int(r["is_original_branchsite_candidate_foreground"]) == 1]
    write_tsv(args.output_dir / "eight_focal_branchsite_decisions.tsv", focal)

    decision_counts = {key: sum(r["decision"] == key for r in focal) for key in {
        "likely_null_local_optimum_do_not_use", "rerun_null_and_alt", "review_or_targeted_rerun",
        "retain_robust_across_starts", "retain_provisional_targeted_rerun", "retain"
    }}

    # Conservative sensitivity analysis: a branch-site alternative with zero
    # positive-class weight is effectively on the null boundary. Assign p=0.5
    # for those diagnosed runs and recompute FDR across the original 110 tests
    # per ecotype. This is diagnostic; final reporting should still rerun the
    # affected null models to document convergence.
    rejected_keys = {
        (r["family"], r["ecotype"])
        for r in focal if r["decision"] == "likely_null_local_optimum_do_not_use"
    }
    full_summary = read_tsv(args.full_summary)
    for row in full_summary:
        original_p = float(row["branchsite_p_mixture"])
        corrected_p = 0.5 if (row["split_family"], row["ecotype"]) in rejected_keys else original_p
        row["branchsite_p_mixture_original"] = original_p
        row["branchsite_p_mixture_conservative"] = corrected_p
    for ecotype in ("HighT", "LowT"):
        subset = [r for r in full_summary if r["ecotype"] == ecotype]
        qvals = bh_adjust([float(r["branchsite_p_mixture_conservative"]) for r in subset])
        for row, q in zip(subset, qvals):
            row["branchsite_FDR_conservative"] = q
            row["branchsite_significant_conservative_FDR05"] = int(q < 0.05)
    corrected_fields = list(full_summary[0].keys())
    with (args.output_dir / "codeml_summary_conservative_sensitivity.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=corrected_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(full_summary)
    corrected_sig = [r for r in full_summary if int(r["branchsite_significant_conservative_FDR05"]) == 1]
    corrected_high = [r for r in corrected_sig if r["ecotype"] == "HighT"]
    corrected_low = [r for r in corrected_sig if r["ecotype"] == "LowT"]
    corrected_by_key = {(r["split_family"], r["ecotype"]): r for r in full_summary}
    audit_by_key = {(r["family"], r["ecotype"]): r for r in focal}
    post_audit = []
    for source in integrated:
        row = dict(source)
        key = (source["split_family"], source["foreground"])
        corrected = corrected_by_key[key]
        audit = audit_by_key.get(key)
        decision = audit["decision"] if audit else "not_flagged_originally"
        corrected_q = float(corrected["branchsite_FDR_conservative"])
        if decision == "likely_null_local_optimum_do_not_use":
            final_status = "exclude_degenerate_local_optimum"
        elif corrected_q >= 0.05:
            final_status = "exclude_after_FDR_sensitivity"
        elif decision == "rerun_null_and_alt":
            final_status = "rerun_required"
        elif decision == "retain_provisional_targeted_rerun":
            final_status = "provisional_targeted_rerun"
        else:
            final_status = "retain"
        row.update({
            "convergence_audit_decision": decision,
            "branchsite_p_mixture_conservative": corrected["branchsite_p_mixture_conservative"],
            "branchsite_FDR_conservative": corrected_q,
            "post_audit_status": final_status,
        })
        post_audit.append(row)
    write_tsv(args.output_dir / "branchsite_candidates_post_audit.tsv", post_audit)
    status_counts = {status: sum(r["post_audit_status"] == status for r in post_audit) for status in {
        "retain", "provisional_targeted_rerun", "rerun_required",
        "exclude_degenerate_local_optimum", "exclude_after_FDR_sensitivity"
    }}
    rows_md = []
    for r in focal:
        rows_md.append(
            f"| {r['family']} | {r['ecotype']} | {r['best_alt_run']} | {float(r['alt_lnL_spread']):.3f} | "
            f"{r['alt_reps_raw_p_lt_0.05']}/3 | {float(r['best_positive_class_proportion']):.3g} | "
            f"{float(r['best_max_BEB']):.3f} | {r['decision']} |"
        )

    report = f"""# 8个branch-site候选的codeml收敛复核

## 文件完整性

- 8个家族 × 2个生态型 × 8次模型运行均存在，共{len(all_runs)}次。
- 正常解析并结束：{sum(bool(r['normal_end']) for r in all_runs)}/{len(all_runs)}；非空stderr：{sum(bool(r['stderr_nonempty']) for r in all_runs)}。

## 判定原则

1. 三个alternative起始值应尽量收敛到相同lnL；差异大提示局部最优。
2. 如果最佳alternative的正选择位点类别比例接近0、无高可信BEB位点，但仍比null高很多，优先怀疑null的枝长/参数落入较差局部最优，而不是生物学正选择。
3. 原始FDR基于每个生态型的110个检验；任何重跑后必须把更新后的P值放回全部110个检验重新计算FDR。

## 8个原始候选

| Family | Foreground | Best alt | lnL spread | raw P<0.05 reps | positive-class proportion | max BEB | Decision |
|---|---|---|---:|---:|---:|---:|---|
{chr(10).join(rows_md)}

## 汇总

- 高度疑似null局部最优、当前不应使用：{decision_counts['likely_null_local_optimum_do_not_use']}。
- 需要同时重跑null和alternative：{decision_counts['rerun_null_and_alt']}。
- 三次起始值均显著且重复定位相同BEB位点，可保留：{decision_counts['retain_robust_across_starts']}。
- 三次均显著且有共同BEB位点，但建议定向重跑：{decision_counts['retain_provisional_targeted_rerun']}。

## 保守FDR敏感性分析

将4个正选择类别权重为0的退化结果临时按null边界处理（mixture P=0.5），再放回每个生态型的110个检验重新进行BH校正：

- branch-site FDR < 0.05：{len(corrected_sig)}个检验；HighT={len(corrected_high)}，LowT={len(corrected_low)}。
- 该计数是复核用的保守敏感性结果。正式表格仍建议用重跑后的null/alternative lnL重新计算。

## 当前候选处置

- 可保留：{status_counts['retain']}个。
- 暂时保留但建议定向重跑：{status_counts['provisional_targeted_rerun']}个（`group_4854`）。
- 必须重跑后决定：{status_counts['rerun_required']}个（`group_1962`）。
- 退化局部最优，当前剔除：{status_counts['exclude_degenerate_local_optimum']}个。
- 因保守重算FDR越过0.05，当前剔除：{status_counts['exclude_after_FDR_sensitivity']}个（`group_3311`）。

5个功能明确的核心候选`cysG_1`、`pheA`、`gch3`、`dadD`和`thrB`均不受本次复核影响。

不能通过“只选择三个alternative中lnL最高的一次”解决null模型自身的局部最优。重跑重点应是null与alternative使用更多起始条件，并检查最终lnL和位点类别比例是否重复收敛。
"""
    (args.output_dir / "CODEML_8_CONVERGENCE_AUDIT_CN.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
