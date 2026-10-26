"""Audit and reproduce a fixed-coefficient coastal AO scaling scenario.

This is a conditional scaling exercise, not a validated process-rate model.
Original abundance models and source tables are never modified. The known
temperature-statistic mismatch is deliberately retained and reported.

Run: python calculate_and_plot.py
Inspect the PNG, then run: python calculate_and_plot.py --export
"""
from pathlib import Path
import hashlib
import json
import platform
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text

import os
if not os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the separate analysis workspace")
_workspace = Path(os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
OUT = _workspace / 'ammonia_oxidation_rate_ESSD/fixed_coefficient_scenario_20261016'
OUT.mkdir(parents=True, exist_ok=True)
ROOT = OUT.parents[1]
AO = OUT.parent
PREV = AO / "rate_abundance_overlay_preview_20261016"
TAG = "nearshore200km_KEEPCHL_SSP585_2080-2100"
STATUS = "exploratory_only_temperature_mismatch_unresolved_not_submission_ready"
UNITS = "nmol N L-1 d-1"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_table(frame, name):
    frame.to_csv(OUT / name, sep="\t", index=False, float_format="%.15g")


def grid_rates(records):
    """Median within study/grid, then median of study medians in each grid."""
    per_study = records.groupby(["grid_id", "study"], as_index=False).agg(
        study_median_rate=("rate", "median"), n_records=("rate", "size"))
    rates = per_study.groupby("grid_id", as_index=False).agg(
        rate_median=("study_median_rate", "median"),
        n_studies=("study", "size"), n_records=("n_records", "sum"))
    return per_study, rates


def prepare():
    original_meta = json.loads((PREV / "summary.json").read_text(encoding="utf-8"))
    if any(Path(p).is_absolute() or ".." in Path(p).parts for p in original_meta["input_hashes"]):
        raise ValueError("input_hashes keys must be workspace-relative paths, with original hashes retained")
    inputs = [PREV / name for name in ["overlay_grid_source_data.tsv",
        "background_abundance_change.tsv", "included_rate_records.tsv", "summary.json"]]
    inputs += [ROOT / p for p in original_meta["input_hashes"]]
    inputs += [AO / "high_temperature_search_20261016/ecotype_temperature_audit.json"]
    inputs = list(dict.fromkeys(inputs))
    before = {str(p.relative_to(ROOT)): digest(p) for p in inputs}
    for p, old_hash in original_meta["input_hashes"].items():
        assert digest(ROOT / p) == old_hash, f"Changed upstream input: {p}"

    g = pd.read_csv(PREV / "background_abundance_change.tsv", sep="\t")
    p = pd.read_csv(PREV / "overlay_grid_source_data.tsv", sep="\t")
    obs = pd.read_csv(PREV / "included_rate_records.tsv", sep="\t")
    audit = pd.read_csv(AO / "coastal_rate_assessment_20261015/raw_selection_audit.tsv", sep="\t")
    selected = audit.loc[audit.valid_coordinates & audit.valid_surface_depth &
        audit.positive_numeric_rate & audit.study.notna() & audit.exact_cell]
    assert set(selected.excel_row) == set(obs.excel_row)
    assert len(g) == 42258 and len(p) == 39 and len(obs) == 100
    assert obs.study.nunique() == 18
    assert g.grid_id.is_unique and p.grid_id.is_unique
    assert not g[["Longitude", "Latitude"]].duplicated().any()
    assert g.grid_depth_m_current.between(-200, 0).all()
    assert g.grid_dist_to_land_km_current.between(0, 200).all()
    assert obs.sample_depth.between(0, 20).all() and obs.rate.gt(0).all()
    assert np.allclose(g.area_weight, np.cos(np.deg2rad(g.Latitude)), atol=1e-14)

    for ec, branch in [("low", "NP_lowT"), ("high", "NP_highT")]:
        base = ROOT / "SSP585" / branch
        training = pd.read_csv(base / "AOA_points_with_BioORACLE_nearshore_filtered.tsv", sep="\t")
        assert training.eps.nunique() == 1
        eps = float(training.eps.iloc[0])
        for period in ["current", "future"]:
            original = pd.read_csv(base / "AOA_future_projection_outputs/GLOBAL_NEARSHORE_MAPS" /
                f"GLOBAL_nearshore_pred_{period}_{TAG}.tsv", sep="\t")
            merged = g.merge(original[["Longitude", "Latitude", f"pred_{period}"]],
                on=["Longitude", "Latitude"], validate="one_to_one")
            assert len(merged) == len(g)
            assert np.allclose(merged[f"{ec}_{period}"],
                10 ** merged[f"pred_{period}"] - eps, rtol=1e-10, atol=1e-12)
    abundance_cols = [f"{e}_{t}" for e in ["low", "high"] for t in ["current", "future"]]
    assert np.isfinite(g[abundance_cols]).all().all()
    assert g[abundance_cols].gt(0).all().all()
    paired = p.merge(g[["grid_id", "Longitude", "Latitude"] + abundance_cols],
        on="grid_id", validate="one_to_one", suffixes=("_obs", "_grid"))
    assert len(paired) == 39
    for col in ["Longitude", "Latitude"] + abundance_cols:
        assert np.allclose(paired[col + "_obs"], paired[col + "_grid"], rtol=1e-12, atol=1e-14)
    per_study, reconstructed = grid_rates(obs)
    comparison = p.merge(reconstructed, on="grid_id", validate="one_to_one", suffixes=("", "_reconstructed"))
    assert np.allclose(comparison.rate_median, comparison.rate_median_reconstructed)
    assert (comparison.n_records == comparison.n_records_reconstructed).all()
    for df in [g, p]:
        df["S_current"] = df.low_current + df.high_current
        df["S_future"] = df.low_future + df.high_future
        df["future_current_factor"] = df.S_future / df.S_current
        df["change_percent"] = 100 * (df.future_current_factor - 1)
    return g, p, obs, per_study, inputs, before


def calculate(g, p, obs, per_study):
    sum_r = float(p.rate_median.sum())
    sum_s = float(p.S_current.sum())
    k = sum_r / sum_s
    p["local_k"] = p.rate_median / p.S_current
    p["scaled_current_rate"] = k * p.S_current
    p["scaled_future_rate"] = k * p.S_future
    p["observed_anchored_future_rate"] = p.rate_median * p.future_current_factor
    p["log10_scaled_over_observed"] = np.log10(p.scaled_current_rate / p.rate_median)
    p["grid_id"] = p.grid_id.astype(int)
    max_idx = p.rate_median.idxmax()
    max_row = p.loc[max_idx]
    rest = p.drop(index=max_idx)
    mean_s0 = float(np.average(g.S_current, weights=g.area_weight))
    mean_s1 = float(np.average(g.S_future, weights=g.area_weight))
    change = 100 * (mean_s1 / mean_s0 - 1)
    assert np.isclose(p.scaled_current_rate.mean(), p.rate_median.mean())
    assert np.isclose(k, np.average(p.local_k, weights=p.S_current))

    definitions = [
        ("sum_ratio_all_grids", "All 39 grids", len(p), k),
        ("sum_ratio_omit_highest_rate", "Omit highest rate", len(rest), float(rest.rate_median.sum() / rest.S_current.sum())),
        ("median_local_ratio", "Median local ratio", len(p), float(p.local_k.median())),
    ]
    scenarios = pd.DataFrame([dict(calibration=key, label=label, n_grids=n, k=coef,
        current_rate=coef*mean_s0, future_rate=coef*mean_s1,
        change_percent=100*(coef*mean_s1/(coef*mean_s0)-1),
        units=UNITS, status=STATUS) for key, label, n, coef in definitions])
    assert np.allclose(scenarios.change_percent, change)

    loo_grid = []
    for row in p.itertuples():
        coef = (sum_r - row.rate_median) / (sum_s - row.S_current)
        loo_grid.append(dict(omitted_grid_id=row.grid_id, k=coef,
            current_rate=coef*mean_s0, future_rate=coef*mean_s1, change_percent=change))
    loo_study = []
    for study in sorted(obs.study.unique()):
        _, remaining = grid_rates(obs.loc[obs.study != study])
        remaining = remaining.merge(p[["grid_id", "S_current"]], on="grid_id", validate="one_to_one")
        coef = float(remaining.rate_median.sum() / remaining.S_current.sum())
        loo_study.append(dict(omitted_study=study, n_grids=len(remaining),
            n_records=int(remaining.n_records.sum()), k=coef,
            current_rate=coef*mean_s0, future_rate=coef*mean_s1, change_percent=change))
    loo_grid = pd.DataFrame(loo_grid)
    loo_study = pd.DataFrame(loo_study)

    g["scaled_current_rate"] = k * g.S_current
    g["scaled_future_rate"] = k * g.S_future
    g["scaled_rate_difference"] = g.scaled_future_rate - g.scaled_current_rate
    assert np.isclose(np.average(g.scaled_current_rate, weights=g.area_weight), k*mean_s0)
    assert np.isclose(np.average(g.scaled_future_rate, weights=g.area_weight), k*mean_s1)
    rho = float(spearmanr(p.S_current, p.rate_median).statistic)
    diag = dict(spearman_rho_descriptive_only=rho,
        in_sample_RMSE_log10=float(np.sqrt(np.mean(p.log10_scaled_over_observed**2))),
        median_absolute_fold_discrepancy=float(10**np.median(abs(p.log10_scaled_over_observed))),
        n_within_factor_2=int((abs(p.log10_scaled_over_observed) <= np.log10(2)).sum()),
        local_k_min=float(p.local_k.min()), local_k_median=float(p.local_k.median()),
        local_k_max=float(p.local_k.max()),
        local_k_max_min_ratio=float(p.local_k.max()/p.local_k.min()),
        maximum_rate_grid_id=int(max_row.grid_id),
        maximum_rate_share_of_sample_sum=float(max_row.rate_median/sum_r),
        maximum_rate_share_of_sample_abundance_sum=float(max_row.S_current/sum_s),
        maximum_rate=float(max_row.rate_median),
        no_holdout_validation_performed=True,
        mean_matching_is_by_construction_not_validation=True)
    meta = dict(status=STATUS, algorithm="ratio-of-sums fixed-coefficient scaling, not OLS",
        n_records=len(obs), n_studies=int(obs.study.nunique()), n_calibration_grids=len(p), n_target_grids=len(g),
        rate_units=UNITS, sum_calibration_rates=sum_r, sum_calibration_S_current=sum_s,
        coefficient_k=k, coefficient_unit="nmol N L-1 d-1 per source relative-abundance unit",
        area_weight_sum=float(g.area_weight.sum()),
        area_weighted_S_current=mean_s0, area_weighted_S_future=mean_s1,
        area_weighted_current_rate=k*mean_s0, area_weighted_future_rate=k*mean_s1,
        relative_change_percent=change, diagnostic=diag,
        alternate_calibrations=scenarios.to_dict("records"),
        leave_one_grid_current_range=[float(loo_grid.current_rate.min()),float(loo_grid.current_rate.max())],
        leave_one_study_current_range=[float(loo_study.current_rate.min()),float(loo_study.current_rate.max())],
        previous_observation_anchored_scenario=dict(current_mean=float(p.rate_median.mean()),
            future_mean=float(p.observed_anchored_future_rate.mean()),
            change_of_means_percent=float(100*(p.observed_anchored_future_rate.mean()/p.rate_median.mean()-1)),
            median_grid_change_percent=float(p.change_percent.median()),
            n_decrease=int((p.change_percent<0).sum()),n_increase=int((p.change_percent>0).sum())),
        uncertainty="No confidence/prediction interval; calibration variants and omission ranges are sensitivity checks only.",
        assumptions=["common scale factor for both ecotypes, all locations, current and future",
            "community AO observations used only to set an empirical reference scale",
            "no independent attribution of observed AO to LowT and HighT",
            "same retained 0.25-degree coastal mask in both periods",
            "cos(latitude) full-cell area weights; no partial-ocean-cell correction",
            "positive quantified 0-20-m observations only, dates pooled, not annual means",
            "no projection, assay, climate, or parameter uncertainty propagated",
            "original abundance temperature mismatch retained at user request"],
        interpretation="The percentage change is algebraically independent of k; it does not test a process-rate response.",
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,
            scipy=scipy.__version__,matplotlib=matplotlib.__version__))
    save_table(p, "source_calibration_39_grids.tsv")
    save_table(g, "source_global_coastal_42258_grids.tsv")
    save_table(obs, "source_included_100_records.tsv")
    save_table(per_study, "source_within_study_grid_medians.tsv")
    save_table(scenarios, "source_global_scenarios.tsv")
    save_table(loo_grid, "sensitivity_leave_one_grid.tsv")
    save_table(loo_study, "sensitivity_leave_one_study.tsv")
    return meta, scenarios


def plot(p, scenarios):
    plt.rcParams.update({"font.family":"Arial", "font.size":7,
        "axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6.5,"ytick.labelsize":6.5,
        "axes.linewidth":1,"axes.spines.top":False,"axes.spines.right":False,
        "legend.frameon":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig = plt.figure(figsize=(180/25.4, 88/25.4))
    a = fig.add_axes([.105,.29,.34,.57])
    b = fig.add_axes([.625,.29,.325,.57])
    # Diagnostic, not held-out validation: both axes in the same volumetric units.
    limits = (.001, 10000)
    a.plot(limits, limits, color=".55", linestyle="--", linewidth=1, zorder=0)
    a.scatter(p.rate_median, p.scaled_current_rate, s=17, color="#3E6D8C",
        edgecolor="white", linewidth=.5, alpha=.85, zorder=2)
    a.set(xscale="log", yscale="log", xlim=limits, ylim=limits)
    ticks = [1e-3,1e-1,1e1,1e3]
    a.set_xticks(ticks)
    a.set_yticks(ticks)
    a.minorticks_off()
    a.set_xlabel("Observed AO rate\n(nmol N L$^{-1}$ d$^{-1}$)", labelpad=4)
    a.set_ylabel("Scaled baseline AO rate\n(nmol N L$^{-1}$ d$^{-1}$)", labelpad=4)
    a.set_title("Baseline check", loc="left", pad=7)
    a.text(.04,.94,"n = 39",transform=a.transAxes,ha="left",va="top",fontsize=7)
    styles = [("#315F7D","o"),("#B0773C","s"),("#7B7B7B","^")]
    handles = []
    for row, (color, marker) in zip(scenarios.itertuples(), styles):
        line, = b.plot([0,1],[row.current_rate,row.future_rate], color=color, marker=marker,
            markersize=4.5, linewidth=1.2, label=row.label)
        handles.append(line)
        for x,y in [(0,row.current_rate),(1,row.future_rate)]:
            b.annotate(f"{y:.1f}",(x,y),xytext=(0,7),textcoords="offset points",ha="center",color=color,fontsize=6.8)
    b.set(xlim=(-.18,1.18),ylim=(0,225))
    b.set_xticks([0,1],["Current","2080–2100"])
    b.set_yticks([0,50,100,150,200])
    b.set_ylabel("Scaled AO rate\n(nmol N L$^{-1}$ d$^{-1}$)",labelpad=5)
    b.set_title("Global coastal scenarios",loc="left",pad=7)
    b.text(.5,.55,"−7.29% by construction",ha="center",transform=b.transAxes,fontsize=7,color=".3")
    for ax, letter in [(a,"a"),(b,"b")]:
        ax.text(-.16,1.07,letter,transform=ax.transAxes,fontweight="bold",fontsize=9)
        ax.tick_params(length=3,width=1)
    legend_obj = fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.53,.055),ncol=3,
        fontsize=6.8,columnspacing=1.6,handlelength=2,handletextpad=.6)
    fig.text(.5,.022,"Fixed-coefficient scenario · provisional inputs",ha="center",fontsize=6.5,color=".4")
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    clipped=[]
    for item in fig.findobj(Text):
        if not item.get_visible() or not item.get_text():
            continue
        box=item.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width+1 or box.y1 > fig.bbox.height+1:
            clipped.append(item.get_text())
    assert not clipped, clipped
    assert not legend_obj.get_window_extent(renderer).overlaps(a.xaxis.label.get_window_extent(renderer))
    stem=OUT / "Supplementary_fixed_coefficient_PREVIEW"
    fig.savefig(stem.with_suffix(".png"),dpi=300)
    if "--export" in sys.argv:
        fig.savefig(stem.with_suffix(".pdf"))
        fig.savefig(stem.with_suffix(".svg"))
        # Raw TIFF avoids the local libtiff LZW encoder's native crash.
        fig.savefig(stem.with_suffix(".tiff"),dpi=600,pil_kwargs={"compression":"raw"})
    plt.close(fig)
    qa=dict(status=STATUS,width_mm=180,height_mm=88,backend="Python/matplotlib",
        font="Arial",minimum_text_points=6.5,text_outside_canvas=clipped,
        panel_a_points=len(p),panel_a_axes="log10; all points included; dashed line is identity",
        panel_a_validation="in-sample diagnostic, not independent validation",
        panel_b_source="source_global_scenarios.tsv",panel_b_y_axis_starts_at_zero=True,
        intervals="none; three deterministic calibration choices are not confidence bounds",
        editable_text="SVG text nodes and PDF TrueType",visual_review="pending")
    (OUT/"figure_qa.json").write_text(json.dumps(qa,indent=2),encoding="utf-8")


def calculation_text(m):
    d=m["diagnostic"]
    lines=[
        "完整计算记录：全球近岸固定换算系数氨氧化情景（探索版）",
        "状态：原丰度预测温度口径尚未修正，不宜直接作为投稿定稿结果。",
        "",
        "1. 数据与筛选",
        "Tang v2 编译速率，保留坐标有效、采样深度0–20 m、正值且可定量、研究标签有效、落入现有近岸网格的记录。",
        "网格定义：距岸≤200 km，海底水深≤200 m；海底水深不是采样深度。不额外使用25-km邻近匹配。",
        "100条记录，18项研究，39个0.25°网格；27条缺少现场温度仍保留。BDL和零值不纳入，不能代表包括低/零速率在内的无偏总体。",
        "同一研究同一网格先取中位数；再在每个网格内对研究中位数取中位数，记为R_i。不是年平均，也非同时测得的两分支组成/速率。",
        "",
        "2. 丰度尺度",
        "从原输出的log10(abundance+epsilon)逆变换为原始相对丰度单位，不能将对数值直接相加。已逐项与原模型输出核对。",
        "S_i,t=L_i,t+H_i,t。L/H为预测相对丰度，不是相互归一化后和为1的分支组成份额。",
        "原始尺度若整体乘100，k会除以100，kS不变；不可混用一份百分比、一份比例。",
        "",
        "3. 原方法：总和之比定标（ratio-of-sums），非最小二乘回归",
        f"sum_i R_i = {m['sum_calibration_rates']:.12f}",
        f"sum_i S_i,0 = {m['sum_calibration_S_current']:.12f}",
        f"k = sum(R_i)/sum(S_i,0) = {m['coefficient_k']:.12f}",
        "k的量纲：nmol N L−1 d−1 / 每原始相对丰度单位。它是经验尺度因子，不是单位细胞生理活性。",
        "每个网格各出现一次，不按原始记录数重复加权；区域较密集的研究仍会有更多网格，不能声称全球代表性或研究完全等权。",
        "k=mean(R)/mean(S)=sum(S_i*k_i)/sum(S_i)，其中k_i=R_i/S_i。它不等于mean(k_i)，也不等于OLS过原点斜率sum(S_i*R_i)/sum(S_i^2)。",
        "这个定义强制39格预测平均值等于39格观测平均值；均值吻合本身不是模型验证。",
        "",
        "4. 全球近岸逐格换算和面积加权",
        "对相同42,258个网格j：Rhat_j,0=k*S_j,0；Rhat_j,1=k*S_j,1。",
        "使用w_j=cos(latitude_j*pi/180)。对同尺寸经纬网格，完整球面网格面积与cos纬度成正比。未校正近岸格中海水面积比例。",
        "Rbar_t=sum(w_j*Rhat_j,t)/sum(w_j)=k*sum(w_j*S_j,t)/sum(w_j)。这不是总通量，也不是全水柱积分。",
        f"sum(w) = {m['area_weight_sum']:.12f}",
        f"Sbar_0 = {m['area_weighted_S_current']:.12f}",
        f"Sbar_1 = {m['area_weighted_S_future']:.12f}",
        f"Rbar_0 = {m['coefficient_k']:.9f} * {m['area_weighted_S_current']:.12f} = {m['area_weighted_current_rate']:.12f} {UNITS}",
        f"Rbar_1 = {m['coefficient_k']:.9f} * {m['area_weighted_S_future']:.12f} = {m['area_weighted_future_rate']:.12f} {UNITS}",
        f"100*(Rbar_1/Rbar_0-1) = {m['relative_change_percent']:.12f}%",
        "",
        "5. 关键恒等式与解释边界",
        "Rbar_1/Rbar_0 = (k*Sbar_1)/(k*Sbar_0) = Sbar_1/Sbar_0。",
        "因此任何空间和时间恒定的正k都会给出相同百分变化。Tang数据只设置纵轴尺度，未独立验证下降方向、幅度或两分支贡献。",
        "采用同一个k相当于对两个分支使用同一系数，并将其他过程贡献/环境效应压缩入尺度；不是分支实测速率之和。",
        "纳入原丰度模型的环境变化仅通过丰度预测进入；未另建单位活性的环境响应，也未显式预测铵底物效应。",
        "",
        "6. 定标诊断和敏感性（不是置信区间）",
        f"最高速率网格={d['maximum_rate_grid_id']}；速率={d['maximum_rate']:.12g}；占39格速率数值和={100*d['maximum_rate_share_of_sample_sum']:.6f}%。",
        f"39格R/S的范围={d['local_k_min']:.12g}至{d['local_k_max']:.12g}；中位数={d['local_k_median']:.12g}。",
        f"观测速率与预测S的描述性Spearman rho={d['spearman_rho_descriptive_only']:.8f}（非验证、无显著性结论）。",
        f"同样本log10 RMSE={d['in_sample_RMSE_log10']:.8f}；绝对倍数偏差中位数={d['median_absolute_fold_discrepancy']:.8f}。",
        f"39格中只有{d['n_within_factor_2']}格的统一系数换算值在观测0.5–2倍内。",
    ]
    for r in m["alternate_calibrations"]:
        lines.append(f"{r['label']}: k={r['k']:.9f}; 当前={r['current_rate']:.9f}; 未来={r['future_rate']:.9f}; 变化={r['change_percent']:.9f}%")
    lines += [
        "遗漏最高速率网格仅用于影响诊断，不表示它是错误记录；主分析未删除它。中位数局地系数仅作不同定标规则敏感性，不替换原主方案。",
        f"逐格遗漏重定标的当前全球均值范围={m['leave_one_grid_current_range']}。",
        f"逐研究遗漏、重新做两级中位数汇总的当前全球均值范围={m['leave_one_study_current_range']}。",
        "这些范围不是预测区间或置信区间，不含原丰度模型、未来环境、实验方法、空间转移和参数误差。",
        "",
        "7. 与先前39格计算的区别",
        "先前：R1_i=Robs_i*S1_i/S0_i，相当于保留各地点不同的k_i；不推广无观测网格。",
        "现在：Rhat0_j=k*S0_j，Rhat1_j=k*S1_j，全域使用一个k；连有观测地点的基线也变成了定标值，不再逐点保留Robs。",
        "两种方法不能把39格变化−30.51%与全球−7.29%只解释为空间范围差异，换算假设也不同。",
        "",
        "8. 论文用途判断",
        "数学和程序正确，可复现条件情景；不能据此认定190.96/177.03为可靠全球平均实测速率，也不能声称全球全部古菌氨氧化下降。",
        "当前版本适合内部诊断；原温度输入问题修正并重算后，若保留，应作为Supplementary Methods及条件情景附图，在Discussion解释，而不是独立功能结果或摘要主结论。",
        "没有修改原丰度模型、原速率表或投稿稿件。输出具体图注、讨论草稿见manuscript_text.txt。",
        "",
        "数据来源：Tang et al. (2023), https://doi.org/10.5194/essd-15-5039-2023；项目使用的v2数据表。",
        "软件、输入哈希和逐格计算见summary.json及source_*.tsv。",
    ]
    (OUT/"完整计算与算法审查.txt").write_text("\n".join(lines)+"\n",encoding="utf-8")


def main():
    g,p,obs,per_study,inputs,before=prepare()
    meta,scenarios=calculate(g,p,obs,per_study)
    meta["input_sha256"]=before
    calculation_text(meta)
    plot(p,scenarios)
    assert before == {str(f.relative_to(ROOT)):digest(f) for f in inputs}
    meta["original_files_unchanged"]=True
    (OUT/"summary.json").write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps({key:meta[key] for key in ["status","coefficient_k","area_weighted_current_rate",
        "area_weighted_future_rate","relative_change_percent","diagnostic","alternate_calibrations",
        "leave_one_study_current_range"]},indent=2,ensure_ascii=False))


if __name__ == "__main__":
    main()
