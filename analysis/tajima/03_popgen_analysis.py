#!/usr/bin/env python
"""
03_popgen_analysis.py
====================================================================
功能：用 scikit-allel 算群体遗传学指标
输入：
  - vcf/HighT_population.vcf.gz
  - vcf/LowT_population.vcf.gz
输出：
  - results/popgen_summary.tsv       全基因组统计
  - results/popgen_windowed.tsv      滑动窗口统计
  - results/per_sample_pi.tsv        每个样本的 π
====================================================================
"""

# Release path adapter: data and outputs live in a separate user workspace.
import os as _release_os
from pathlib import Path as _ReleasePath
if not _release_os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the analysis workspace")
_release_root = _ReleasePath(_release_os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
_release_os.chdir(_release_root / 'tajima')


import os
import sys
import numpy as np
import pandas as pd
import allel
from pathlib import Path

# -------------------- 配置 --------------------
HIGHT_VCF = "vcf/HighT_population.vcf.gz"
LOWT_VCF  = "vcf/LowT_population.vcf.gz"
WINDOW_SIZE = 10000   # 滑动窗口 10kb
WINDOW_STEP = 5000    # 步长 5kb
MIN_ALT_FRAC = 0.05   # 低于这个频率的 allele 不计

Path("results").mkdir(exist_ok=True)

# -------------------- 工具函数 --------------------
def load_population_vcf(vcf_path, label):
    """读入 VCF，返回 allele counts array 和元数据"""
    print(f"\n[{label}] 读入 VCF: {vcf_path}")
    callset = allel.read_vcf(
        vcf_path,
        fields=['variants/CHROM', 'variants/POS', 'variants/REF', 'variants/ALT',
                'calldata/GT', 'calldata/AD', 'samples']
    )
    if callset is None:
        raise ValueError(f"无法读取 {vcf_path}")
    
    samples = callset['samples']
    pos = callset['variants/POS']
    chrom = callset['variants/CHROM']
    gt = allel.GenotypeArray(callset['calldata/GT'])
    
    print(f"  samples: {len(samples)}")
    print(f"  SNVs:    {len(pos)}")
    print(f"  chroms:  {len(np.unique(chrom))}")
    
    # Allele counts (2D: sites × alleles)
    ac = gt.count_alleles()
    
    # 筛选双等位 + 分离位点（segregating sites）
    is_biallelic = ac.is_biallelic_01()
    is_segregating = ac.is_segregating()
    mask = is_biallelic & is_segregating
    
    pos_f = pos[mask]
    chrom_f = chrom[mask]
    ac_f = ac.compress(mask, axis=0)
    gt_f = gt.compress(mask, axis=0)
    
    print(f"  Biallelic segregating SNVs: {len(pos_f)}")
    
    return {
        'samples': samples,
        'chrom': chrom_f,
        'pos': pos_f,
        'ac': ac_f,
        'gt': gt_f,
    }


def compute_global_stats(data, label):
    """算全基因组 π, θ, Tajima's D, ..."""
    ac = data['ac']
    pos = data['pos']
    n_samples = len(data['samples'])
    
    # Sequence diversity (π)
    pi = allel.sequence_diversity(pos, ac)
    # Watterson's theta
    theta = allel.watterson_theta(pos, ac)
    # Tajima's D
    D = allel.tajima_d(ac)
    
    result = {
        'clade': label,
        'n_samples': n_samples,
        'n_SNVs': len(pos),
        'pi': pi,
        'theta_W': theta,
        'Tajima_D': D,
    }
    
    print(f"\n[{label}] Global statistics:")
    for k, v in result.items():
        if isinstance(v, float):
            print(f"  {k:15s} = {v:.6f}")
        else:
            print(f"  {k:15s} = {v}")
    
    return result


def compute_windowed_stats(data, label, window_size, window_step):
    """滑动窗口 π + Tajima's D（按 contig 分别算）"""
    records = []
    pos = data['pos']
    ac = data['ac']
    chrom = data['chrom']
    
    # 按 contig 分
    for contig in np.unique(chrom):
        mask = chrom == contig
        if mask.sum() < 10:  # 窗口太少
            continue
        p_c = pos[mask]
        ac_c = ac.compress(mask, axis=0)
        
        max_pos = p_c.max()
        windows = np.array([
            [s, min(s + window_size - 1, max_pos)]
            for s in range(1, max_pos, window_step)
        ])
        
        # π per window
        pi_w, _, _, _ = allel.windowed_diversity(
            p_c, ac_c, windows=windows
        )
        # Tajima's D per window
        D_w, _, _ = allel.windowed_tajima_d(
            p_c, ac_c, windows=windows
        )
        
        for i, w in enumerate(windows):
            records.append({
                'clade': label,
                'chrom': contig,
                'start': w[0],
                'end': w[1],
                'pi': pi_w[i],
                'Tajima_D': D_w[i],
            })
    
    return pd.DataFrame(records)


def compute_per_sample_pi(data, label):
    """给每个样本算它的 π（用它的 ALT allele frequency）"""
    gt = data['gt']
    pos = data['pos']
    samples = data['samples']
    
    records = []
    for i, s in enumerate(samples):
        # 对单个样本：π 约等于该样本作为 pool 时的期望异质度
        # 这里用一个简化估计：2 * f * (1-f) 的均值，f 是该样本的 ALT 频率
        # （对 pooled-continuous VCF 更合适，但如果 GT 是 0/0, 0/1, 1/1 也能近似）
        
        # 取单样本 allele counts
        single_gt = gt[:, i:i+1, :]
        single_ac = allel.GenotypeArray(single_gt).count_alleles()
        
        # 该样本的 π
        pi_s = allel.sequence_diversity(pos, single_ac)
        
        records.append({
            'sample': s,
            'clade': label,
            'pi': pi_s,
        })
    
    return pd.DataFrame(records)


# -------------------- 主流程 --------------------
print("=" * 60)
print(" 群体遗传学分析 (scikit-allel)")
print("=" * 60)

# 读入两个 clade
highT = load_population_vcf(HIGHT_VCF, "HighT")
lowT  = load_population_vcf(LOWT_VCF,  "LowT")

# -------------------- 全基因组统计 --------------------
global_stats = pd.DataFrame([
    compute_global_stats(highT, "HighT"),
    compute_global_stats(lowT,  "LowT"),
])
global_stats.to_csv("results/popgen_summary.tsv", sep="\t", index=False)

# -------------------- 滑动窗口统计 --------------------
print("\n" + "=" * 60)
print(" 滑动窗口分析")
print("=" * 60)

windowed = pd.concat([
    compute_windowed_stats(highT, "HighT", WINDOW_SIZE, WINDOW_STEP),
    compute_windowed_stats(lowT,  "LowT",  WINDOW_SIZE, WINDOW_STEP),
], ignore_index=True)
windowed.to_csv("results/popgen_windowed.tsv", sep="\t", index=False)
print(f"  窗口数: {len(windowed)} (HighT: {(windowed.clade=='HighT').sum()}, LowT: {(windowed.clade=='LowT').sum()})")

# -------------------- Per-sample π --------------------
print("\n" + "=" * 60)
print(" Per-sample π")
print("=" * 60)

per_sample = pd.concat([
    compute_per_sample_pi(highT, "HighT"),
    compute_per_sample_pi(lowT,  "LowT"),
], ignore_index=True)
per_sample.to_csv("results/per_sample_pi.tsv", sep="\t", index=False)
print(per_sample.groupby('clade')['pi'].describe())

# -------------------- 统计检验 --------------------
from scipy.stats import mannwhitneyu, ttest_1samp

print("\n" + "=" * 60)
print(" 统计检验")
print("=" * 60)

# 1. π (HighT vs LowT) - Mann-Whitney U on windowed values
highT_pi_win = windowed[(windowed.clade=="HighT") & windowed.pi.notna()].pi
lowT_pi_win  = windowed[(windowed.clade=="LowT")  & windowed.pi.notna()].pi
u, p = mannwhitneyu(lowT_pi_win, highT_pi_win, alternative='greater')
print(f"\nπ: LowT > HighT ?")
print(f"  Mann-Whitney U = {u:.0f}, p = {p:.2e}")
print(f"  LowT median π = {lowT_pi_win.median():.5f}")
print(f"  HighT median π = {highT_pi_win.median():.5f}")

# 2. Tajima's D (HighT vs LowT)
highT_D_win = windowed[(windowed.clade=="HighT") & windowed.Tajima_D.notna()].Tajima_D
lowT_D_win  = windowed[(windowed.clade=="LowT")  & windowed.Tajima_D.notna()].Tajima_D
u, p = mannwhitneyu(highT_D_win, lowT_D_win, alternative='less')
print(f"\nTajima's D: HighT < LowT ?")
print(f"  Mann-Whitney U = {u:.0f}, p = {p:.2e}")
print(f"  HighT median D = {highT_D_win.median():.3f}")
print(f"  LowT median D = {lowT_D_win.median():.3f}")

# 3. HighT 的 D 是否显著 < 0?
t, p = ttest_1samp(highT_D_win, 0, alternative='less')
print(f"\nHighT Tajima's D < 0 ?")
print(f"  one-sample t-test: t = {t:.2f}, p = {p:.2e}")
print(f"  (显著负 → 支持近期正选择/扩张假说)")

# 4. LowT 的 D 是否显著不等于 0?
t, p = ttest_1samp(lowT_D_win, 0)
print(f"\nLowT Tajima's D ≠ 0 ?")
print(f"  one-sample t-test: t = {t:.2f}, p = {p:.2e}")
print(f"  (不显著 → 支持稳态中性演化)")

# -------------------- 保存统计检验 --------------------
with open("results/statistical_tests.txt", "w") as f:
    f.write("Population genetics statistical tests\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"HighT median π (windowed) = {highT_pi_win.median():.5f}\n")
    f.write(f"LowT  median π (windowed) = {lowT_pi_win.median():.5f}\n")
    f.write(f"HighT median Tajima's D (windowed) = {highT_D_win.median():.3f}\n")
    f.write(f"LowT  median Tajima's D (windowed) = {lowT_D_win.median():.3f}\n")

print("\n" + "=" * 60)
print(" 完成！结果保存在 results/")
print("=" * 60)
print("  results/popgen_summary.tsv       - 全基因组统计")
print("  results/popgen_windowed.tsv      - 滑动窗口")
print("  results/per_sample_pi.tsv        - per-sample π")
print("  results/statistical_tests.txt    - 统计检验汇总")
