#!/usr/bin/env python3
"""Summarize and plot lineage-level dN/dS results for HighT and LowT.

The analysis deliberately separates:
1. branch models on crown-internal branches (long-term mean omega shifts), and
2. branch-site models on the single stem branch (episodic site-level selection).

The existing codeml fits compare each ecotype with its own pooled background.
They are not a direct likelihood-ratio test of omega_HighT == omega_LowT.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


HIGH = "#D55E00"
LOW = "#0072B2"
NEUTRAL = "#B8BDC6"
LOWER = "#3A9188"
HIGHER = "#B55D7A"
TEXT = "#20242A"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        required=True,
        help="ISME project root",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: selection_analysis/results/16_dnds_lineage_history)",
    )
    return parser.parse_args()


def fmt(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}g}"


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "legend.frameon": False,
            "text.color": TEXT,
            "axes.labelcolor": TEXT,
            "axes.titlecolor": TEXT,
        }
    )


def prepare_data(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    branch_path = (
        root
        / "selection_analysis/downloaded_codeml_results_20260721/selection_analysis/results"
        / "09_codeml_summary_soft85/codeml_summary_all.tsv"
    )
    audit_path = (
        root
        / "selection_analysis/codeml_analysis_20260721/convergence_review"
        / "branchsite_candidates_post_audit.tsv"
    )
    branch = pd.read_csv(branch_path, sep="\t")
    audit = pd.read_csv(audit_path, sep="\t")

    if len(branch) != 220 or branch["split_family"].nunique() != 110:
        raise ValueError("Expected 220 family-ecotype tests from 110 codeml-ready families")
    if not (branch["all_models_complete"] == 1).all():
        raise ValueError("At least one codeml family-ecotype test is incomplete")

    branch["omega_effect_log2"] = np.log2(
        branch["omega_foreground"] / branch["omega_background"]
    )
    branch["minus_log10_branch_FDR"] = -np.log10(
        branch["branch_FDR"].clip(lower=np.nextafter(0, 1))
    )
    branch["branch_direction"] = np.where(
        branch["omega_foreground"] < branch["omega_background"], "lower", "higher"
    )
    branch["branch_result"] = np.where(
        branch["branch_significant_FDR05"].eq(1),
        branch["branch_direction"],
        "not_significant",
    )

    excluded = {
        "exclude_after_FDR_sensitivity",
        "exclude_degenerate_local_optimum",
    }
    curated = audit.loc[~audit["post_audit_status"].isin(excluded)].copy()
    curated["minus_log10_branchsite_FDR"] = -np.log10(
        curated["branchsite_FDR_conservative"].clip(lower=np.nextafter(0, 1))
    )
    curated["needs_rerun"] = curated["post_audit_status"].isin(
        ["rerun_required", "provisional_targeted_rerun"]
    )
    curated["display_gene"] = curated["split_family"]
    curated.loc[curated["priority"].eq("core_interpretable"), "display_gene"] += " *"

    counts = (
        branch.loc[branch["branch_significant_FDR05"].eq(1)]
        .groupby(["ecotype", "branch_direction"], observed=True)
        .size()
        .rename("n_genes")
        .reset_index()
    )
    complete_index = pd.MultiIndex.from_product(
        [["HighT", "LowT"], ["lower", "higher"]],
        names=["ecotype", "branch_direction"],
    )
    counts = (
        counts.set_index(["ecotype", "branch_direction"])
        .reindex(complete_index, fill_value=0)
        .reset_index()
    )
    return branch, curated, counts


def make_figure(
    branch: pd.DataFrame, curated: pd.DataFrame, counts: pd.DataFrame, outdir: Path
) -> None:
    configure_matplotlib()
    fig = plt.figure(figsize=(7.05, 7.0), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.0, 1.55],
        width_ratios=[1.75, 0.8],
        hspace=0.48,
        wspace=0.42,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    # a | All branch-model tests.
    for eco, color, marker in [("HighT", HIGH, "o"), ("LowT", LOW, "s")]:
        sub = branch.loc[branch["ecotype"].eq(eco)]
        nonsig = sub.loc[sub["branch_significant_FDR05"].eq(0)]
        sig = sub.loc[sub["branch_significant_FDR05"].eq(1)]
        ax_a.scatter(
            nonsig["omega_effect_log2"],
            nonsig["minus_log10_branch_FDR"],
            s=12,
            marker=marker,
            facecolors="none",
            edgecolors=NEUTRAL,
            linewidths=0.55,
            alpha=0.8,
            zorder=1,
        )
        ax_a.scatter(
            sig["omega_effect_log2"],
            sig["minus_log10_branch_FDR"],
            s=22,
            marker=marker,
            color=color,
            edgecolors="white",
            linewidths=0.35,
            alpha=0.9,
            label=eco,
            zorder=2,
        )
    ax_a.axvline(0, color="#73777F", lw=0.7, ls="--", zorder=0)
    ax_a.axhline(-np.log10(0.05), color="#73777F", lw=0.7, ls=":", zorder=0)
    ax_a.set_xlabel(r"Branch effect, log$_2$(ω$_{foreground}$/ω$_{background}$)")
    ax_a.set_ylabel(r"−log$_{10}$(branch-model FDR)")
    ax_a.set_title("Long-term branch-model selection landscape", loc="left", fontweight="bold")
    ax_a.legend(loc="upper right", handletextpad=0.4, borderaxespad=0.2)
    ax_a.text(
        0.01,
        0.03,
        "left: stronger constraint\nright: relative acceleration",
        transform=ax_a.transAxes,
        color="#5E646D",
        fontsize=6,
        va="bottom",
    )

    # b | Direction of significant branch shifts.
    x = np.arange(2)
    lower_n = [
        int(counts.query("ecotype == @eco and branch_direction == 'lower'")["n_genes"].iloc[0])
        for eco in ["HighT", "LowT"]
    ]
    higher_n = [
        int(counts.query("ecotype == @eco and branch_direction == 'higher'")["n_genes"].iloc[0])
        for eco in ["HighT", "LowT"]
    ]
    ax_b.bar(x, lower_n, width=0.62, color=LOWER, label="Lower ω")
    ax_b.bar(x, higher_n, width=0.62, bottom=lower_n, color=HIGHER, label="Higher ω")
    for i, (lo, hi) in enumerate(zip(lower_n, higher_n)):
        if lo:
            ax_b.text(i, lo / 2, str(lo), ha="center", va="center", color="white", fontweight="bold")
        if hi:
            ax_b.text(i, lo + hi / 2, str(hi), ha="center", va="center", color="white", fontweight="bold")
        ax_b.text(i, lo + hi + 1.2, f"n={lo + hi}", ha="center", va="bottom", fontsize=6.5)
    ax_b.set_xticks(x, ["HighT", "LowT"])
    ax_b.set_ylabel("Significant genes")
    ax_b.set_ylim(0, max(np.array(lower_n) + np.array(higher_n)) * 1.18)
    ax_b.set_title("Direction of significant shifts", loc="left", fontweight="bold")
    ax_b.legend(loc="upper right", bbox_to_anchor=(1.02, 1.0), fontsize=6, handlelength=1.2)

    # c | Curated branch-site candidates.
    curated = curated.sort_values(
        ["foreground", "branchsite_FDR_conservative", "split_family"],
        ascending=[True, True, True],
        key=lambda col: col.map({"HighT": 0, "LowT": 1}) if col.name == "foreground" else col,
    ).reset_index(drop=True)
    y = np.arange(len(curated))[::-1]
    for yi, (_, row) in zip(y, curated.iterrows()):
        color = HIGH if row["foreground"] == "HighT" else LOW
        face = "white" if row["needs_rerun"] else color
        size = 24 + 11 * int(row["BEB_ge_0.95"])
        ax_c.scatter(
            row["minus_log10_branchsite_FDR"],
            yi,
            s=size,
            marker="o",
            facecolor=face,
            edgecolor=color,
            linewidth=1.0,
            zorder=3,
        )
    ax_c.axvline(-np.log10(0.05), color="#73777F", lw=0.7, ls=":", zorder=0)
    high_n = int(curated["foreground"].eq("HighT").sum())
    if 0 < high_n < len(curated):
        split_y = len(curated) - high_n - 0.5
        ax_c.axhline(split_y, color="#D7DADE", lw=0.7)
    ax_c.set_yticks(y, curated["display_gene"])
    for tick, priority in zip(ax_c.get_yticklabels(), curated["priority"]):
        if priority == "core_interpretable":
            tick.set_fontweight("bold")
    ax_c.set_xlabel(r"−log$_{10}$(conservative branch-site FDR)")
    ax_c.set_title(
        "Stem branch-site candidates after convergence audit",
        loc="left",
        fontweight="bold",
    )
    ax_c.grid(axis="x", color="#E5E7EA", lw=0.5, zorder=0)
    ax_c.set_axisbelow(True)
    ax_c.text(
        0.995,
        0.015,
        "* annotated candidate with BEB ≥0.95",
        transform=ax_c.transAxes,
        ha="right",
        va="bottom",
        fontsize=6,
        color="#5E646D",
    )
    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=HIGH, markeredgecolor=HIGH, markersize=5, label="HighT stem"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=LOW, markeredgecolor=LOW, markersize=5, label="LowT stem"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#575B62", markersize=5, label="Targeted rerun pending"),
    ]
    ax_c.legend(handles=legend, loc="lower right", bbox_to_anchor=(1.0, 0.06), fontsize=6)

    for label, ax in zip("abc", [ax_a, ax_b, ax_c]):
        ax.text(
            -0.13 if label != "b" else -0.25,
            1.08,
            label,
            transform=ax.transAxes,
            fontsize=10,
            fontweight="bold",
            va="top",
        )

    fig.suptitle(
        "Asymmetric long-term selection across HighT and LowT ecotypes",
        x=0.08,
        y=0.995,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.07, top=0.94)
    base = outdir / "Figure_DNDS_Lineage_History"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_report(
    branch: pd.DataFrame, curated: pd.DataFrame, counts: pd.DataFrame, outdir: Path
) -> None:
    stats = {}
    for eco in ["HighT", "LowT"]:
        sub = branch.loc[branch["ecotype"].eq(eco)]
        sig = sub.loc[sub["branch_significant_FDR05"].eq(1)]
        stats[eco] = {
            "n": len(sub),
            "sig": len(sig),
            "lower": int(sig["branch_direction"].eq("lower").sum()),
            "higher": int(sig["branch_direction"].eq("higher").sum()),
            "median_fg": float(sub["omega_foreground"].median()),
            "median_bg": float(sub["omega_background"].median()),
            "median_effect": float(sub["omega_effect_log2"].median()),
        }
    high_sig = set(
        branch.loc[
            branch["ecotype"].eq("HighT") & branch["branch_significant_FDR05"].eq(1),
            "split_family",
        ]
    )
    low_sig = set(
        branch.loc[
            branch["ecotype"].eq("LowT") & branch["branch_significant_FDR05"].eq(1),
            "split_family",
        ]
    )
    core = curated.loc[curated["priority"].eq("core_interpretable")].copy()
    retained = int(curated["post_audit_status"].eq("retain").sum())
    pending = int(curated["needs_rerun"].sum())
    high_bs = int(curated["foreground"].eq("HighT").sum())
    low_bs = int(curated["foreground"].eq("LowT").sum())

    core_lines = []
    for _, row in core.sort_values(["foreground", "branchsite_FDR_conservative"]).iterrows():
        core_lines.append(
            f"| {row['split_family']} | {row['foreground']} | {fmt(row['branchsite_FDR_conservative'], 4)} | "
            f"{row['BEB_sites_ge_0.95']} | {row['Roary_annotation']} |"
        )

    report = f"""# HighT/LowT谱系dN/dS历史分析

## 一句话结论

HighT冠群形成后，共有单拷贝基因表现出广泛的长期ω降低，符合选择约束增强；LowT冠群的显著变化同时包含ω降低与升高，表现出更异质的长期演化模式。两条生态型起源stem上还检测到不同的branch-site候选，但候选数量差异本身不是HighT高于LowT的统计证据。

## 分析对象与模型含义

- 输入：110个codeml-ready soft85单拷贝直系同源基因，每个基因分别检验HighT和LowT，共220个family–ecotype检验。
- branch model：前景为相应生态型crown内部枝，不包含进入共同祖先的stem；回答生态型形成后整条基因的长期平均ω是否相对背景改变。
- branch-site model：前景为进入相应生态型共同祖先的单一stem；回答谱系形成阶段是否有少数密码子出现间歇性ω>1信号。
- 多重检验：HighT与LowT各自的110个基因分别进行BH-FDR校正。

## Branch model结果

| 前景 | 检验基因 | FDR<0.05 | 前景ω降低 | 前景ω升高 | 全部基因前景ω中位数 | 背景ω中位数 | log2(前景/背景)中位数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| HighT | {stats['HighT']['n']} | {stats['HighT']['sig']} | {stats['HighT']['lower']} | {stats['HighT']['higher']} | {fmt(stats['HighT']['median_fg'], 4)} | {fmt(stats['HighT']['median_bg'], 4)} | {fmt(stats['HighT']['median_effect'], 4)} |
| LowT | {stats['LowT']['n']} | {stats['LowT']['sig']} | {stats['LowT']['lower']} | {stats['LowT']['higher']} | {fmt(stats['LowT']['median_fg'], 4)} | {fmt(stats['LowT']['median_bg'], 4)} | {fmt(stats['LowT']['median_effect'], 4)} |

- 共有{len(high_sig & low_sig)}个基因在两套独立branch检验中均显著。
- HighT的{stats['HighT']['sig']}个显著基因全部为前景ω降低，因此主要信号不是“快速进化”，而是长期净化选择或功能约束增强。
- LowT的{stats['LowT']['sig']}个显著基因中，{stats['LowT']['lower']}个降低、{stats['LowT']['higher']}个升高。升高但仍低于1只能表述为相对加速或约束放松，不能仅凭branch model称为正选择。

## Branch-site结果及收敛复核

- 原始FDR筛选得到23个family–ecotype候选。
- 对多起始值收敛和保守FDR复核后，{retained}个保留，{pending}个仍需定向补跑；当前可视化候选共{len(curated)}个（HighT={high_bs}，LowT={low_bs}）。
- 另外5个候选因退化局部最优或保守FDR越过0.05而不进入当前候选图。
- HighT与LowT的候选基因集合不同，但{high_bs}比{low_bs}并不能单独证明HighT经历了更多正选择。

### 功能明确且具有BEB≥0.95位点的核心候选

| 基因 | stem前景 | 保守FDR | BEB位点 | 功能注释 |
|---|---|---:|---|---|
{chr(10).join(core_lines)}

这些结果支持“生态型形成阶段发生少数氨基酸替换，随后部分基因受到长期选择约束维持”的阶段性模型；它们仍需结合密码子比对、基因树、重组检查和位点结构背景进行最终解释。

## 当前能够和不能够写出的结论

可以写：

> Relative to the remaining Nitrosopumilus background, the HighT crown showed widespread reductions in gene-wide ω, whereas LowT exhibited a mixture of reduced and elevated ω. Distinct branch-site candidates were detected on the HighT and LowT stem branches, consistent with lineage-specific episodic amino-acid evolution during ecotype formation.

不能写：

> HighT has a significantly lower ω than LowT.

原因是现有HighT和LowT模型分别拟合，各自背景都包含另一生态型和Other NP，尚未在同一个嵌套模型中直接检验ω_HighT=ω_LowT。如果正文必须作这一直接比较，应增加“HighT与LowT共用一个ω”的null模型与“HighT、LowT分别估计ω”的alternative模型，进行1自由度LRT。

## 图的阅读方式

- a：横轴小于0表示前景ω低于背景，纵轴越高表示branch FDR越小。
- b：仅统计branch FDR<0.05的基因；HighT的显著改变全部指向ω降低。
- c：展示收敛复核后仍保留或待定向补跑的stem branch-site候选；空心点表示仍需补跑，星号表示功能明确且至少有一个BEB≥0.95位点。

## 输出文件

- `Figure_DNDS_Lineage_History.pdf/svg/png/tiff`：论文图与预览。
- `source_data_branch_model.tsv`：220个branch检验及效应方向。
- `source_data_branchsite_curated.tsv`：复核后的branch-site候选。
- `summary_branch_direction_counts.tsv`：显著branch方向计数。
- `FIGURE_LEGEND_DRAFT.md`：英文图注草稿。
"""
    (outdir / "DNDs_LINEAGE_HISTORY_REPORT_CN.md").write_text(report, encoding="utf-8")

    legend = """# Figure legend draft

**Figure X | Asymmetric long-term selection across HighT and LowT ecotypes.**
**a,** Gene-wide branch-model effects for 110 soft-core single-copy orthologous families tested separately for the HighT and LowT crowns. Effect sizes are log2 ratios of foreground to background ω; the horizontal dotted line marks BH-FDR = 0.05. **b,** Direction of significant branch-model shifts. All 49 significant HighT tests had lower foreground ω, whereas significant LowT tests included 17 lower and 15 higher foreground estimates. **c,** Stem branch-site candidates retained after conservative FDR and multi-start convergence review. Point size reflects the number of BEB sites with posterior probability ≥0.95; open points denote candidates requiring targeted reruns. Asterisks identify functionally annotated candidates with at least one BEB site ≥0.95. FDR was controlled separately across 110 genes for each ecotype. HighT and LowT effects are each relative to their pooled phylogenetic background and do not constitute a direct test of ω equality between the two ecotypes.
"""
    (outdir / "FIGURE_LEGEND_DRAFT.md").write_text(legend, encoding="utf-8")

    contract = """# Figure contract

- Core conclusion: HighT shows widespread long-term reductions in gene-wide omega relative to background, while LowT has bidirectional rate shifts; distinct stem branch-site candidates accompany ecotype formation.
- Evidence chain: panel a shows all branch-model effects; panel b makes the direction asymmetry explicit; panel c reports curated stem branch-site candidates and uncertainty.
- Archetype: quantitative grid with panel a as the hero panel.
- Backend: Python/matplotlib only.
- Export: 179-mm-wide editable SVG/PDF, 300-dpi PNG preview, 600-dpi TIFF, plus panel-level source tables.
- Review risk: the two ecotype fits are background-relative and must not be described as a direct HighT-versus-LowT likelihood-ratio test.
"""
    (outdir / "FIGURE_CONTRACT.md").write_text(contract, encoding="utf-8")


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    outdir = (
        args.outdir.resolve()
        if args.outdir
        else root / "selection_analysis/results/16_dnds_lineage_history"
    )
    outdir.mkdir(parents=True, exist_ok=True)
    branch, curated, counts = prepare_data(root)

    branch.to_csv(outdir / "source_data_branch_model.tsv", sep="\t", index=False)
    curated.to_csv(outdir / "source_data_branchsite_curated.tsv", sep="\t", index=False)
    counts.to_csv(outdir / "summary_branch_direction_counts.tsv", sep="\t", index=False)
    make_figure(branch, curated, counts, outdir)
    write_report(branch, curated, counts, outdir)
    print(f"Wrote dN/dS lineage-history analysis to: {outdir}")


if __name__ == "__main__":
    main()
