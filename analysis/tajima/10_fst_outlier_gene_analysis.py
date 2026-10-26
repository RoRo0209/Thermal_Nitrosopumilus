#!/usr/bin/env python
"""
10_fst_outlier_gene_analysis.py
====================================================================
Fst outlier + 基因注释 + 功能分类分析

输入:
  - results/fst_windowed.tsv        (per-window Fst, 来自 07_fst_analysis.py)
  - refs/LowT_representative.gff    (reference 基因组注释，GFF3)
  - LowT_kofam_clean.tsv            (parse_kofam_output.py 的输出)
  - (可选) candidate_genes.tsv      (你前面筛的候选基因列表)

输出:
  - results/fst_outlier_genes_top5.tsv   top 5% outlier 里的基因
  - results/fst_outlier_genes_top1.tsv   top 1% outlier 里的基因
  - results/fst_candidate_overlap.tsv    候选基因在 outlier 里的分布
  - figures/fst_distribution.pdf         Fst 分布 + 阈值线
  - figures/fst_manhattan.pdf             沿基因组的 Fst 分布
  - figures/fst_outlier_functions.pdf    Outlier gene 的功能分类

用法:
  python 10_fst_outlier_gene_analysis.py
====================================================================
"""

# Release path adapter: data and outputs live in a separate user workspace.
import os as _release_os
from pathlib import Path as _ReleasePath
if not _release_os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the analysis workspace")
_release_root = _ReleasePath(_release_os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
_release_os.chdir(_release_root / 'tajima')


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import re
import argparse

Path("results").mkdir(exist_ok=True)
Path("figures").mkdir(exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--fst', default='results/fst_windowed.tsv')
parser.add_argument('--gff', default='refs/LowT_representative.gff')
parser.add_argument('--kofam', default='LowT_kofam_clean.tsv')
parser.add_argument('--candidates', default='candidate_genes.tsv',
                    help='可选：候选基因列表')
parser.add_argument('--top_percent_main', type=float, default=5.0,
                    help='主分析阈值 (top %%)')
parser.add_argument('--top_percent_strict', type=float, default=1.0,
                    help='严格阈值 (top %%)')
args = parser.parse_args()

# ==========================================================
# Step 1: 读 Fst 数据
# ==========================================================
print("=" * 60)
print(" Step 1: 读 per-window Fst")
print("=" * 60)

fst = pd.read_csv(args.fst, sep='\t')
fst = fst.dropna(subset=['Fst']).copy()
# 确保 Fst >= 0（某些估计量会给负值，视作 0）
fst['Fst_raw'] = fst['Fst']
fst['Fst'] = fst['Fst'].clip(lower=0)

print(f"  Windows: {len(fst)}")
print(f"  Contigs: {fst['chrom'].nunique()}")
print(f"  Fst 分布:")
print(f"    mean   = {fst['Fst'].mean():.3f}")
print(f"    median = {fst['Fst'].median():.3f}")
print(f"    min    = {fst['Fst'].min():.3f}")
print(f"    max    = {fst['Fst'].max():.3f}")

# 定义 outlier 阈值
thr_main   = fst['Fst'].quantile(1 - args.top_percent_main / 100)
thr_strict = fst['Fst'].quantile(1 - args.top_percent_strict / 100)

print(f"\n  Top {args.top_percent_main}% threshold:   Fst >= {thr_main:.3f}")
print(f"  Top {args.top_percent_strict}% threshold:   Fst >= {thr_strict:.3f}")

fst['is_outlier_main']   = fst['Fst'] >= thr_main
fst['is_outlier_strict'] = fst['Fst'] >= thr_strict

print(f"  N outlier (top {args.top_percent_main}%):   {fst['is_outlier_main'].sum()}")
print(f"  N outlier (top {args.top_percent_strict}%): {fst['is_outlier_strict'].sum()}")


# ==========================================================
# Step 2: 解析 GFF
# ==========================================================
print("\n" + "=" * 60)
print(" Step 2: 解析 GFF")
print("=" * 60)

def parse_gff_attrs(attr_str):
    """Parse GFF3 attribute column"""
    d = {}
    for kv in attr_str.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k.strip()] = v.strip()
    return d

gff_rows = []
with open(args.gff) as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            continue
        if parts[2] != 'CDS':
            continue
        attrs = parse_gff_attrs(parts[8])
        gene_id = attrs.get('ID', attrs.get('locus_tag', ''))
        # Prokka 通常用 ID=PREFIX_XXXXX
        gff_rows.append({
            'chrom': parts[0],
            'start': int(parts[3]),
            'end':   int(parts[4]),
            'strand': parts[6],
            'gene_id': gene_id,
            'product': attrs.get('product', ''),
            'locus_tag': attrs.get('locus_tag', ''),
        })

gff_df = pd.DataFrame(gff_rows)
print(f"  CDS 数: {len(gff_df)}")
print(f"  Contigs in GFF: {gff_df['chrom'].nunique()}")
print(f"  GFF 里的 gene_id 样例: {gff_df['gene_id'].iloc[0]}")


# ==========================================================
# Step 3: 读 KEGG/KO 注释
# ==========================================================
print("\n" + "=" * 60)
print(" Step 3: 读 KO 注释")
print("=" * 60)

kofam = pd.read_csv(args.kofam, sep='\t')
print(f"  KO 表 rows: {len(kofam)}")
print(f"  KO 表列: {list(kofam.columns)}")
print(f"  kofam gene_id 样例: {kofam['gene_id'].iloc[0]}")

# 诊断 gene_id 一致性
gff_ids = set(gff_df['gene_id'])
kofam_ids = set(kofam['gene_id'])
overlap = gff_ids & kofam_ids
print(f"\n  GFF gene_id 数:     {len(gff_ids)}")
print(f"  kofam gene_id 数:   {len(kofam_ids)}")
print(f"  交集:               {len(overlap)}")

if len(overlap) < 0.5 * min(len(gff_ids), len(kofam_ids)):
    print("  ⚠ 交集太小！gene_id 可能不匹配")
    print(f"    GFF 样例:   {list(gff_ids)[:3]}")
    print(f"    kofam 样例: {list(kofam_ids)[:3]}")
    # 尝试自动修复：去掉 kofam 前缀
    print("\n  尝试自动去前缀...")
    # 猜测前缀：找 kofam 和 GFF 的最长公共后缀
    sample_kofam = list(kofam_ids)[0]
    sample_gff = list(gff_ids)[0]
    # 如果 gff_id 是 kofam_id 的后缀
    if sample_kofam.endswith(sample_gff[-min(20, len(sample_gff)):]):
        # 找到分界点
        for i in range(len(sample_kofam)):
            suffix = sample_kofam[i:]
            if suffix in gff_ids:
                prefix_to_strip = sample_kofam[:i]
                print(f"    识别前缀: '{prefix_to_strip}'")
                kofam['gene_id'] = kofam['gene_id'].str.replace(
                    f'^{re.escape(prefix_to_strip)}', '', regex=True)
                break
        # 重新检查
        kofam_ids2 = set(kofam['gene_id'])
        overlap2 = gff_ids & kofam_ids2
        print(f"    修复后交集: {len(overlap2)}")


# ==========================================================
# Step 4: Overlap outlier windows with genes
# ==========================================================
print("\n" + "=" * 60)
print(" Step 4: 把 outlier windows 映射到基因")
print("=" * 60)

def windows_to_genes(windows_df, gff_df):
    """返回: window 和 gene 的 overlap 关系"""
    records = []
    for contig, win_sub in windows_df.groupby('chrom'):
        gene_sub = gff_df[gff_df['chrom'] == contig]
        if gene_sub.empty:
            continue
        for _, w in win_sub.iterrows():
            # 和 window 重叠的基因：gene_end >= window_start AND gene_start <= window_end
            hits = gene_sub[
                (gene_sub['end']   >= w['start']) &
                (gene_sub['start'] <= w['end'])
            ]
            for _, g in hits.iterrows():
                records.append({
                    'chrom': w['chrom'],
                    'window_start': w['start'],
                    'window_end': w['end'],
                    'Fst': w['Fst'],
                    'gene_id': g['gene_id'],
                    'gene_start': g['start'],
                    'gene_end': g['end'],
                    'product': g['product'],
                })
    return pd.DataFrame(records)

out_main_genes = windows_to_genes(fst[fst['is_outlier_main']], gff_df)
out_strict_genes = windows_to_genes(fst[fst['is_outlier_strict']], gff_df)

print(f"  Top 5% windows overlap genes: {len(out_main_genes)} records, "
      f"{out_main_genes['gene_id'].nunique()} unique genes")
print(f"  Top 1% windows overlap genes: {len(out_strict_genes)} records, "
      f"{out_strict_genes['gene_id'].nunique()} unique genes")


# ==========================================================
# Step 5: Merge KO annotation
# ==========================================================
print("\n" + "=" * 60)
print(" Step 5: Merge KO 注释")
print("=" * 60)

def merge_kofam(gene_df, kofam_df):
    m = gene_df.merge(
        kofam_df[['gene_id', 'KO', 'description']],
        on='gene_id', how='left'
    )
    # 如果 kofam 没有，就用 GFF 的 product
    m['annotation'] = m['description'].fillna(m['product'])
    return m

out_main = merge_kofam(out_main_genes, kofam)
out_strict = merge_kofam(out_strict_genes, kofam)

# 每个 gene 只保留它所在窗口中 Fst 最大值的那一条（去重）
def dedup_gene(df):
    return (df.sort_values('Fst', ascending=False)
              .drop_duplicates('gene_id', keep='first'))

out_main_unique = dedup_gene(out_main)
out_strict_unique = dedup_gene(out_strict)

# 保存
out_main_unique.to_csv(f'results/fst_outlier_genes_top{int(args.top_percent_main)}.tsv',
                       sep='\t', index=False)
out_strict_unique.to_csv(f'results/fst_outlier_genes_top{int(args.top_percent_strict)}.tsv',
                         sep='\t', index=False)

print(f"  Top {args.top_percent_main}% unique genes: {len(out_main_unique)}")
print(f"  Top {args.top_percent_strict}% unique genes: {len(out_strict_unique)}")


# ==========================================================
# Step 6: Top outlier genes 展示
# ==========================================================
print("\n" + "=" * 60)
print(f" Step 6: Top 30 outlier genes (top {args.top_percent_main}%)")
print("=" * 60)

top30 = out_main_unique.head(30)
for _, r in top30.iterrows():
    ko = r['KO'] if pd.notna(r['KO']) else '-'
    ann = r['annotation'] if pd.notna(r['annotation']) else 'hypothetical'
    print(f"  Fst={r['Fst']:.3f}  {r['gene_id']:30s}  {ko:8s}  {ann[:60]}")


# ==========================================================
# Step 7: 候选基因闭环（如果提供）
# ==========================================================
if Path(args.candidates).exists():
    print("\n" + "=" * 60)
    print(" Step 7: 候选基因闭环验证")
    print("=" * 60)
    
    cand = pd.read_csv(args.candidates, sep='\t')
    print(f"  候选基因数: {len(cand)}")
    
    # 按 symbol 或 KO 匹配 outlier
    if 'symbol' in cand.columns:
        out_main_unique['symbol_guess'] = out_main_unique['description'].str.extract(
            r'^([a-zA-Z][\w\d\-]+);', expand=False)
        hit_by_symbol = out_main_unique[
            out_main_unique['symbol_guess'].isin(cand['symbol'])
        ]
        print(f"  Symbol 匹配 outlier: {len(hit_by_symbol)}")
    
    if 'KO' in cand.columns:
        hit_by_ko = out_main_unique[out_main_unique['KO'].isin(cand['KO'])]
        print(f"  KO 匹配 outlier: {len(hit_by_ko)}")
        if len(hit_by_ko) > 0:
            print("\n  候选基因在 outlier 里的条目:")
            for _, r in hit_by_ko.iterrows():
                print(f"    Fst={r['Fst']:.3f}  {r['KO']}  {r['annotation'][:60]}")
    
    # 保存
    if 'KO' in cand.columns:
        overlap_df = hit_by_ko.copy()
        overlap_df.to_csv('results/fst_candidate_overlap.tsv', sep='\t', index=False)
else:
    print(f"\n  候选基因文件 {args.candidates} 不存在，跳过闭环验证")


# ==========================================================
# Step 8: 出图
# ==========================================================
print("\n" + "=" * 60)
print(" Step 8: 出图")
print("=" * 60)

# --- Figure 1: Fst distribution ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(fst['Fst'], bins=50, color='#B0C4DE', edgecolor='black', alpha=0.8)
ax.axvline(thr_main, color='#E8A07A', linestyle='--', linewidth=2,
           label=f'Top {args.top_percent_main}% threshold = {thr_main:.3f}')
ax.axvline(thr_strict, color='#D47746', linestyle='--', linewidth=2,
           label=f'Top {args.top_percent_strict}% threshold = {thr_strict:.3f}')
ax.axvline(fst['Fst'].median(), color='grey', linestyle=':',
           label=f'Median = {fst["Fst"].median():.3f}')
ax.set_xlabel('Per-window Fst')
ax.set_ylabel('Number of windows')
ax.set_title(f'Distribution of genome-wide Fst (n={len(fst)} windows)')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('figures/fst_distribution.pdf', bbox_inches='tight', dpi=300)
plt.savefig('figures/fst_distribution.png', bbox_inches='tight', dpi=150)
print("  figures/fst_distribution.pdf")

# --- Figure 2: Manhattan-like plot ---
# 累计坐标
fst_sorted = fst.sort_values(['chrom', 'start']).copy()
chrom_order = fst_sorted['chrom'].unique()
chrom_sizes = fst_sorted.groupby('chrom')['end'].max().to_dict()

cumsum = 0
chrom_offsets = {}
for c in chrom_order:
    chrom_offsets[c] = cumsum
    cumsum += chrom_sizes[c]

fst_sorted['x_cum'] = fst_sorted.apply(
    lambda r: chrom_offsets[r['chrom']] + (r['start'] + r['end']) / 2, axis=1)

fig, ax = plt.subplots(figsize=(14, 5))

# 按 contig 着色（交替）
colors = ['#4A8FA8', '#9CCFDE']
for i, c in enumerate(chrom_order):
    sub = fst_sorted[fst_sorted['chrom'] == c]
    ax.scatter(sub['x_cum'], sub['Fst'],
               color=colors[i % 2], s=10, alpha=0.6, edgecolor='none')

# 标 outlier
outliers = fst_sorted[fst_sorted['is_outlier_main']]
ax.scatter(outliers['x_cum'], outliers['Fst'],
           color='#E8A07A', s=30, alpha=0.9, edgecolor='black', linewidth=0.3,
           label=f'Top {args.top_percent_main}% outliers (n={len(outliers)})')

# 阈值线
ax.axhline(thr_main, color='#D47746', linestyle='--', linewidth=1,
           label=f'Threshold = {thr_main:.3f}')

ax.set_xlabel('Cumulative genome position')
ax.set_ylabel('Per-window Fst')
ax.set_title('Genome-wide Fst landscape (Manhattan-like)')
ax.legend(loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('figures/fst_manhattan.pdf', bbox_inches='tight', dpi=300)
plt.savefig('figures/fst_manhattan.png', bbox_inches='tight', dpi=150)
print("  figures/fst_manhattan.pdf")

# --- Figure 3: Outlier functional classification ---
def classify_function(desc):
    if pd.isna(desc):
        return 'Unknown / hypothetical'
    d = desc.lower()
    if any(k in d for k in ['monooxygenase', 'amoa', 'amob', 'amoc',
                            'nitrogen', 'nitr', 'urease', 'urea',
                            'ammoni', 'glutam', 'glnb', 'glnk']):
        return 'Nitrogen metabolism'
    if any(k in d for k in ['transport', 'permease', 'abc', 'channel',
                            'symporter', 'antiporter']):
        return 'Transport'
    if any(k in d for k in ['chaperon', 'heat shock', 'hsp', 'stress',
                            'groel', 'dnak', 'protease']):
        return 'Stress/chaperone'
    if any(k in d for k in ['ribos', 'trna', 'rrna', 'translat',
                            'aminoacyl', 'elongat', 'initiation']):
        return 'Translation'
    if any(k in d for k in ['dna', 'replicat', 'repair', 'helicase',
                            'topoisomerase', 'polymerase ii', 'polymerase iii']):
        return 'DNA repair/replication'
    if any(k in d for k in ['lipid', 'fatty', 'membrane', 'phospho lip',
                            'isoprenoid']):
        return 'Membrane/lipid'
    if any(k in d for k in ['dehydrogenase', 'oxidoreductase', 'reductase',
                            'oxidase', 'metabolism', 'kinase', 'synthase',
                            'transferase', 'hydrolase']):
        return 'Central metabolism'
    if any(k in d for k in ['transcription', 'regulator', 'sigma',
                            'two-component', 'response regul']):
        return 'Regulation'
    return 'Other annotated'

out_main_unique['category'] = out_main_unique['annotation'].apply(classify_function)

cat_counts = out_main_unique['category'].value_counts()

fig, ax = plt.subplots(figsize=(8, 5))
cat_counts.sort_values().plot(kind='barh', ax=ax,
                               color='#4A8FA8', edgecolor='black')
ax.set_xlabel('Number of outlier genes')
ax.set_title(f'Functional categories of top {args.top_percent_main}% Fst outlier genes')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('figures/fst_outlier_functions.pdf', bbox_inches='tight', dpi=300)
plt.savefig('figures/fst_outlier_functions.png', bbox_inches='tight', dpi=150)
print("  figures/fst_outlier_functions.pdf")


# ==========================================================
# Summary
# ==========================================================
print("\n" + "=" * 60)
print(" 完成")
print("=" * 60)
print(f"Total windows analyzed:      {len(fst)}")
print(f"Outlier windows (top {args.top_percent_main}%):   {fst['is_outlier_main'].sum()}")
print(f"Outlier windows (top {args.top_percent_strict}%): {fst['is_outlier_strict'].sum()}")
print(f"Outlier genes (top {args.top_percent_main}%):     {out_main_unique['gene_id'].nunique()}")
print(f"Outlier genes (top {args.top_percent_strict}%):   {out_strict_unique['gene_id'].nunique()}")
print(f"\n输出文件:")
print(f"  results/fst_outlier_genes_top{int(args.top_percent_main)}.tsv")
print(f"  results/fst_outlier_genes_top{int(args.top_percent_strict)}.tsv")
print(f"  figures/fst_distribution.pdf")
print(f"  figures/fst_manhattan.pdf")
print(f"  figures/fst_outlier_functions.pdf")
