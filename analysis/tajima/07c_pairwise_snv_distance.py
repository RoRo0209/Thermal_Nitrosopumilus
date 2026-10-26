#!/usr/bin/env python
"""
07c_pairwise_snv_distance.py
====================================================================
适用于 pooled-continuous VCF 的 pairwise 距离计算

不用 Fst（对 1/1 格式不适用），改用:
  - SNV distance: sample i 和 j consensus 不同的位点比例
  - 这是经典的微生物群体遗传距离度量（相当于 1 - ANI）
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

highT_samples = [s.strip() for s in open("HighT_samples_passing.txt") if s.strip()]
lowT_samples  = [s.strip() for s in open("LowT_samples_passing.txt")  if s.strip()]

print("读入 joint VCF...")
callset = allel.read_vcf(JOINT_VCF,
                         fields=['variants/POS', 'calldata/GT', 'samples'])
samples = list(callset['samples'])
highT_idx = [i for i, s in enumerate(samples) if s in highT_samples]
lowT_idx  = [i for i, s in enumerate(samples) if s in lowT_samples]

gt = allel.GenotypeArray(callset['calldata/GT'])
pos = callset['variants/POS']

# 取出每个样本在每个位点的 "consensus allele" (对 haploid / pooled-continuous 模式)
# 0 = ref, 1 = alt, -1 = missing
# gt shape: (n_sites, n_samples, ploidy)
# 对于 haploid (ploidy=1) 或 pooled (都是 0/0 或 1/1)，取第一个即可
gt_array = gt.values  # 3D: sites × samples × ploidy

# 只保留至少一个样本是 variant 的位点
is_variant = gt.count_alleles().is_variant()
pos_v = pos[is_variant]
gt_v = gt_array[is_variant]  # shape: (n_var_sites, n_samples, ploidy)

# 取每个样本的 consensus allele（第一个非 -1 的值）
# 对 haploid/pooled 情况，ploidy 维度上两个值应该一样
# 用 max across ploidy axis（-1 会被当作 -1，但 0/1 会正确）
consensus = gt_v[:, :, 0]  # 取第一个 allele call
# 处理 missing：用 -1 表示
# 如果 allele=-1 说明该位点该样本 missing

print(f"\nVariant sites: {len(pos_v)}")
print(f"Consensus matrix shape: {consensus.shape}")

# -------------------- 计算 pairwise SNV distance --------------------
print("\n计算 pairwise SNV distance...")

n_samples = len(samples)
dist_matrix = np.full((n_samples, n_samples), np.nan)
n_sites_used = np.zeros((n_samples, n_samples), dtype=int)

MIN_SHARED_SITES = 100

for i in range(n_samples):
    for j in range(i, n_samples):
        if i == j:
            dist_matrix[i, j] = 0
            continue
        
        ci = consensus[:, i]
        cj = consensus[:, j]
        
        # 只用两个样本都有 call 的位点
        both_called = (ci >= 0) & (cj >= 0)
        n_shared = both_called.sum()
        
        if n_shared < MIN_SHARED_SITES:
            continue
        
        # 不同的位点数
        n_diff = (ci[both_called] != cj[both_called]).sum()
        
        # SNV distance = 不同比例
        dist = n_diff / n_shared
        
        dist_matrix[i, j] = dist_matrix[j, i] = dist
        n_sites_used[i, j] = n_sites_used[j, i] = n_shared

# -------------------- 摘要 --------------------
highT_pairs = [(i, j) for i in highT_idx for j in highT_idx if i < j]
lowT_pairs  = [(i, j) for i in lowT_idx  for j in lowT_idx  if i < j]
between     = [(i, j) for i in highT_idx for j in lowT_idx]

def get_valid_dists(pair_list):
    return [dist_matrix[i, j] for i, j in pair_list if not np.isnan(dist_matrix[i, j])]

highT_d = get_valid_dists(highT_pairs)
lowT_d  = get_valid_dists(lowT_pairs)
between_d = get_valid_dists(between)

print("\n" + "="*60)
print(" Pairwise SNV distance (proportion of differing sites)")
print("="*60)
print(f"Within HighT ({len(highT_d)} pairs):")
print(f"  median = {np.median(highT_d):.4f}  ({np.median(highT_d)*100:.2f}%)")
print(f"  mean   = {np.mean(highT_d):.4f}")
print(f"  range  = [{min(highT_d):.4f}, {max(highT_d):.4f}]")

print(f"\nWithin LowT ({len(lowT_d)} pairs):")
print(f"  median = {np.median(lowT_d):.4f}  ({np.median(lowT_d)*100:.2f}%)")
print(f"  mean   = {np.mean(lowT_d):.4f}")
print(f"  range  = [{min(lowT_d):.4f}, {max(lowT_d):.4f}]")

print(f"\nBetween HighT-LowT ({len(between_d)} pairs):")
print(f"  median = {np.median(between_d):.4f}  ({np.median(between_d)*100:.2f}%)")
print(f"  mean   = {np.mean(between_d):.4f}")
print(f"  range  = [{min(between_d):.4f}, {max(between_d):.4f}]")

# Fold differences
print(f"\n对比 (between / within):")
print(f"  Between / Within HighT: {np.median(between_d)/np.median(highT_d):.1f}×")
print(f"  Between / Within LowT:  {np.median(between_d)/np.median(lowT_d):.1f}×")

# -------------------- 保存 --------------------
dist_df = pd.DataFrame(dist_matrix, index=samples, columns=samples)
dist_df.to_csv("results/pairwise_snv_distance.tsv", sep="\t")

# -------------------- 统计检验 --------------------
from scipy.stats import mannwhitneyu
u1, p1 = mannwhitneyu(between_d, highT_d, alternative='greater')
u2, p2 = mannwhitneyu(between_d, lowT_d, alternative='greater')
print(f"\n统计检验 (Mann-Whitney U):")
print(f"  Between > Within HighT: p = {p1:.2e}")
print(f"  Between > Within LowT:  p = {p2:.2e}")

# -------------------- 出图 --------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel 1: Heatmap
ordered_samples = [s for s in samples if s in highT_samples] + \
                  [s for s in samples if s in lowT_samples]
dist_df_ord = dist_df.loc[ordered_samples, ordered_samples]

ax = axes[0]
plot_mat = dist_df_ord.values.copy()
np.fill_diagonal(plot_mat, np.nan)

vmax_use = max(0.02, np.nanpercentile(plot_mat, 98))
sns.heatmap(plot_mat, cmap='RdYlBu_r', vmin=0, vmax=vmax_use,
            square=True, cbar_kws={'label': 'SNV distance'},
            ax=ax, xticklabels=False, yticklabels=False)

n_highT_actual = sum(1 for s in ordered_samples if s in highT_samples)
ax.axhline(n_highT_actual, color='black', linewidth=2.5)
ax.axvline(n_highT_actual, color='black', linewidth=2.5)
ax.text(n_highT_actual/2, -2, 'HighT', ha='center',
        color='#D47746', fontweight='bold', fontsize=11)
ax.text(n_highT_actual + (len(ordered_samples)-n_highT_actual)/2, -2, 'LowT',
        ha='center', color='#4A8FA8', fontweight='bold', fontsize=11)
ax.set_title('Pairwise consensus differences among variant sites')

# Panel 2: Distribution
ax = axes[1]
data_plot = [highT_d, lowT_d, between_d]
labels_plot = [f'Within HighT\nn={len(highT_d)}',
               f'Within LowT\nn={len(lowT_d)}',
               f'HighT vs LowT\nn={len(between_d)}']
colors_plot = ['#E8A07A', '#9CCFDE', '#C4B8CE']

parts = ax.violinplot(data_plot, showmedians=True, showextrema=False, widths=0.7)
for pc, c in zip(parts['bodies'], colors_plot):
    pc.set_facecolor(c); pc.set_alpha(0.8); pc.set_edgecolor('black')

# Jittered scatter
for i, d in enumerate(data_plot):
    x_jitter = np.random.normal(i+1, 0.05, len(d))
    ax.scatter(x_jitter, d, s=10, alpha=0.4, color='black', zorder=3)

ax.set_xticks([1, 2, 3])
ax.set_xticklabels(labels_plot)
ax.set_ylabel('SNV distance (proportion of differing sites)')
ax.set_title('Within vs between-ecotype distance')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 标注统计
y_max = ax.get_ylim()[1]
ax.text(3.3, y_max * 0.9,
        f'Between > HighT: p = {p1:.1e}\nBetween > LowT: p = {p2:.1e}',
        fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='grey'))

plt.tight_layout()
plt.savefig('figures/snv_distance_analysis.pdf', bbox_inches='tight', dpi=300)
plt.savefig('figures/snv_distance_analysis.png', bbox_inches='tight', dpi=200)
print("\n图保存: figures/snv_distance_analysis.{pdf,png}")

# -------------------- Summary --------------------
with open('results/snv_distance_summary.txt', 'w') as f:
    f.write("Pairwise SNV Distance Summary\n")
    f.write("="*60 + "\n\n")
    f.write(f"Within HighT ({len(highT_d)} pairs):\n")
    f.write(f"  median = {np.median(highT_d):.4f} ({np.median(highT_d)*100:.2f}%)\n")
    f.write(f"Within LowT ({len(lowT_d)} pairs):\n")
    f.write(f"  median = {np.median(lowT_d):.4f} ({np.median(lowT_d)*100:.2f}%)\n")
    f.write(f"Between HighT-LowT ({len(between_d)} pairs):\n")
    f.write(f"  median = {np.median(between_d):.4f} ({np.median(between_d)*100:.2f}%)\n")
    f.write(f"\nMann-Whitney U tests:\n")
    f.write(f"  Between > Within HighT: p = {p1:.2e}\n")
    f.write(f"  Between > Within LowT:  p = {p2:.2e}\n")

print("\nSummary: results/snv_distance_summary.txt")
