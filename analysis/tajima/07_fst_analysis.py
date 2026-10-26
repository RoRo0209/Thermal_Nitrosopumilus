#!/usr/bin/env python
"""
07_fst_analysis.py
====================================================================
Fst 分析：
  1. Genome-wide HighT vs LowT Fst (Hudson's estimator)
  2. Per-window Fst
  3. Pairwise Fst matrix (每个 sample pair)
====================================================================
"""

# Release path adapter: data and outputs live in a separate user workspace.
import os as _release_os
from pathlib import Path as _ReleasePath
if not _release_os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the analysis workspace")
_release_root = _ReleasePath(_release_os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
_release_os.chdir(_release_root / 'tajima')


import allel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

Path("results").mkdir(exist_ok=True)
Path("figures").mkdir(exist_ok=True)

JOINT_VCF = "vcf_shared/HighT_LowT_joint.vcf.gz"

# 读 populations 映射表
highT_samples = [s.strip() for s in open("HighT_samples_passing.txt") if s.strip()]
lowT_samples  = [s.strip() for s in open("LowT_samples_passing.txt")  if s.strip()]

print(f"HighT samples: {len(highT_samples)}")
print(f"LowT samples:  {len(lowT_samples)}")

# -------------------- 读 joint VCF --------------------
print("\n读入 joint VCF...")
callset = allel.read_vcf(JOINT_VCF,
                         fields=['variants/POS', 'variants/CHROM',
                                 'calldata/GT', 'samples'])
samples = list(callset['samples'])
print(f"  joint samples: {len(samples)}")

# 找出 HighT 和 LowT 在 joint VCF 中的索引
highT_idx = [i for i, s in enumerate(samples) if s in highT_samples]
lowT_idx  = [i for i, s in enumerate(samples) if s in lowT_samples]

print(f"  HighT in joint: {len(highT_idx)}")
print(f"  LowT in joint:  {len(lowT_idx)}")

# -------------------- 计算 allele counts 每个 pop --------------------
gt = allel.GenotypeArray(callset['calldata/GT'])
pos = callset['variants/POS']

ac_H = gt.count_alleles(subpop=highT_idx)
ac_L = gt.count_alleles(subpop=lowT_idx)

# 筛选双等位 + 在两 pop 中都至少有一个非缺失
is_biallelic = ac_H.is_biallelic() & ac_L.is_biallelic()
# 至少总和有变异
is_seg = (ac_H.sum(axis=1) > 0) & (ac_L.sum(axis=1) > 0)
mask = is_biallelic & is_seg

pos_f = pos[mask]
ac_H_f = ac_H.compress(mask, axis=0)
ac_L_f = ac_L.compress(mask, axis=0)
print(f"  biallelic + covered sites: {mask.sum()}")

# -------------------- Genome-wide Fst (Hudson) --------------------
num, den = allel.hudson_fst(ac_H_f, ac_L_f)
fst_global = np.nansum(num) / np.nansum(den)
print(f"\n===== Genome-wide Hudson's Fst =====")
print(f"  Fst (HighT vs LowT) = {fst_global:.3f}")

# No generic species-boundary interpretation is assigned to this exploratory statistic.

# -------------------- Per-window Fst --------------------
print("\n计算 per-window Fst (10kb 窗口)...")
# 按 contig 分别算
chrom = callset['variants/CHROM'][mask]
window_records = []
for c in np.unique(chrom):
    cmask = chrom == c
    if cmask.sum() < 10:
        continue
    p_c = pos_f[cmask]
    ac_H_c = ac_H_f.compress(cmask, axis=0)
    ac_L_c = ac_L_f.compress(cmask, axis=0)
    max_pos = p_c.max()
    windows = np.array([[s, min(s+10000-1, max_pos)]
                        for s in range(1, max_pos, 5000)])
    # Official scikit-allel provides Hudson components and windowed_statistic.
    # Aggregate the components before dividing, preserving the Hudson estimator.
    num_c, den_c = allel.hudson_fst(ac_H_c, ac_L_c)
    num_w, _, _ = allel.windowed_statistic(p_c, num_c, statistic=np.nansum, windows=windows)
    den_w, _, _ = allel.windowed_statistic(p_c, den_c, statistic=np.nansum, windows=windows)
    fst_w = np.divide(num_w, den_w, out=np.full_like(num_w, np.nan, dtype=float), where=den_w > 0)
    for i, w in enumerate(windows):
        window_records.append({
            'chrom': c, 'start': w[0], 'end': w[1], 'Fst': fst_w[i]
        })

windowed_fst = pd.DataFrame(window_records)
windowed_fst.to_csv("results/fst_windowed.tsv", sep="\t", index=False)
print(f"  窗口数: {len(windowed_fst)}")
print(f"  Fst median: {windowed_fst.Fst.median():.3f}")
print(f"  Fst 95% range: [{windowed_fst.Fst.quantile(0.025):.3f}, {windowed_fst.Fst.quantile(0.975):.3f}]")

# -------------------- Pairwise Fst matrix --------------------
print("\n计算 pairwise Fst matrix...")
n_total = len(samples)
fst_matrix = np.zeros((n_total, n_total))

for i in range(n_total):
    for j in range(i+1, n_total):
        ac_i = gt.count_alleles(subpop=[i])
        ac_j = gt.count_alleles(subpop=[j])
        
        m = ac_i.is_biallelic() & ac_j.is_biallelic() & \
            (ac_i.sum(axis=1) > 0) & (ac_j.sum(axis=1) > 0)
        if m.sum() < 50:
            fst_matrix[i, j] = fst_matrix[j, i] = np.nan
            continue
        
        num, den = allel.hudson_fst(ac_i.compress(m, axis=0), ac_j.compress(m, axis=0))
        fst_pair = np.nansum(num) / np.nansum(den) if np.nansum(den) > 0 else np.nan
        fst_matrix[i, j] = fst_matrix[j, i] = fst_pair

# 转成 DataFrame
fst_df = pd.DataFrame(fst_matrix, index=samples, columns=samples)
fst_df.to_csv("results/pairwise_fst_matrix.tsv", sep="\t")

# 打印一些 pairwise Fst 摘要
print("\n===== Pairwise Fst 摘要 =====")
# HighT 内部
highT_pairs = [fst_matrix[i,j] for i in highT_idx for j in highT_idx if i<j]
# LowT 内部
lowT_pairs = [fst_matrix[i,j] for i in lowT_idx for j in lowT_idx if i<j]
# Between
between = [fst_matrix[i,j] for i in highT_idx for j in lowT_idx]

print(f"HighT 内部 (n pairs = {len(highT_pairs)}):")
print(f"  mean = {np.nanmean(highT_pairs):.3f}, median = {np.nanmedian(highT_pairs):.3f}")
print(f"LowT 内部 (n pairs = {len(lowT_pairs)}):")
print(f"  mean = {np.nanmean(lowT_pairs):.3f}, median = {np.nanmedian(lowT_pairs):.3f}")
print(f"HighT vs LowT (n pairs = {len(between)}):")
print(f"  mean = {np.nanmean(between):.3f}, median = {np.nanmedian(between):.3f}")

# -------------------- 出图 --------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Panel 1: Pairwise Fst heatmap (按 clade 排序)
ordered_samples = [s for s in samples if s in highT_samples] + \
                  [s for s in samples if s in lowT_samples]
fst_df_ord = fst_df.loc[ordered_samples, ordered_samples]
np.fill_diagonal(fst_df_ord.values, np.nan)

ax = axes[0]
sns.heatmap(fst_df_ord, cmap='RdYlBu_r', vmin=0, vmax=min(1, np.nanmax(fst_matrix)),
            square=True, cbar_kws={'label': 'Fst'}, ax=ax,
            xticklabels=False, yticklabels=False)
# 画出 HighT/LowT 分界线
n_highT_actual = sum(1 for s in ordered_samples if s in highT_samples)
ax.axhline(n_highT_actual, color='black', linewidth=2)
ax.axvline(n_highT_actual, color='black', linewidth=2)
ax.text(n_highT_actual/2, -1, 'HighT', ha='center', color='#E8A07A', fontweight='bold')
ax.text(n_highT_actual + (len(ordered_samples)-n_highT_actual)/2, -1, 'LowT',
        ha='center', color='#9CCFDE', fontweight='bold')
ax.set_title(f'Pairwise Fst (genome-wide Fst = {fst_global:.3f})')

# Panel 2: Fst distributions
ax = axes[1]
data_plot = [highT_pairs, lowT_pairs, between]
labels_plot = [f'Within\nHighT\nn={len(highT_pairs)}',
               f'Within\nLowT\nn={len(lowT_pairs)}',
               f'HighT\nvs LowT\nn={len(between)}']
colors_plot = ['#E8A07A', '#9CCFDE', '#CCCCCC']

finite_groups = [np.asarray(x)[np.isfinite(x)] for x in data_plot]
valid_positions = [i + 1 for i, values in enumerate(finite_groups) if len(values) > 1 and np.ptp(values) > 0]
if valid_positions:
    parts = ax.violinplot([finite_groups[i - 1] for i in valid_positions], positions=valid_positions,
                         showmedians=True, showextrema=False)
    for pc, i in zip(parts['bodies'], valid_positions):
        pc.set_facecolor(colors_plot[i - 1]); pc.set_alpha(0.8); pc.set_edgecolor('black')
for i, values in enumerate(finite_groups, 1):
    if i not in valid_positions:
        if len(values):
            ax.scatter([i], [np.median(values)], color=colors_plot[i - 1])
        else:
            ax.text(i, 0, 'Not estimable', ha='center', fontsize=8)

ax.set_xticks([1, 2, 3])
ax.set_xticklabels(labels_plot)
ax.set_ylabel('Pairwise Fst')
ax.axhline(fst_global, color='red', linestyle='--', label=f'Genome Fst = {fst_global:.3f}')
ax.legend()
ax.set_title('Fst distribution (within vs between ecotypes)')

plt.tight_layout()
plt.savefig('figures/fst_analysis.pdf', bbox_inches='tight', dpi=300)
plt.savefig('figures/fst_analysis.png', bbox_inches='tight', dpi=200)
print("\n图保存: figures/fst_analysis.{pdf,png}")

# -------------------- 写入 summary --------------------
with open('results/fst_summary.txt', 'w') as f:
    f.write("Fst Analysis Summary\n")
    f.write("="*60 + "\n\n")
    f.write(f"Genome-wide Hudson's Fst (HighT vs LowT): {fst_global:.3f}\n")
    f.write(f"Interpretation: {interp}\n\n")
    f.write(f"Pairwise Fst:\n")
    f.write(f"  Within HighT:  mean = {np.nanmean(highT_pairs):.3f}, "
            f"median = {np.nanmedian(highT_pairs):.3f}\n")
    f.write(f"  Within LowT:   mean = {np.nanmean(lowT_pairs):.3f}, "
            f"median = {np.nanmedian(lowT_pairs):.3f}\n")
    f.write(f"  Between HighT-LowT: mean = {np.nanmean(between):.3f}, "
            f"median = {np.nanmedian(between):.3f}\n\n")
    f.write(f"Per-window Fst (10kb):\n")
    f.write(f"  median = {windowed_fst.Fst.median():.3f}\n")
    f.write(f"  95% range = [{windowed_fst.Fst.quantile(0.025):.3f}, "
            f"{windowed_fst.Fst.quantile(0.975):.3f}]\n")
print("Summary 保存: results/fst_summary.txt")
