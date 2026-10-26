#!/bin/bash
: "${THERMAL_PROJECT_ROOT:?Set THERMAL_PROJECT_ROOT to the analysis workspace}"
cd "${THERMAL_PROJECT_ROOT}/abundance_wca_np"

###############################################################################
# OGA amoA Ecotype 分析 v3 —— 修正NP_other + 去重
#
# 修正内容:
#   1. NP = HighT + LowT + NP_other (之前漏了NP_other)
#   2. 去除重复样本 (如TARA_100_DCM出现4次)
#   3. NP_other纳入known ecotypes和total计算
#   4. 图例和图表包含NP_other
#
# 用法: bash run_oga_analysis_v3.sh /path/to/ref_all_amoA.fasta
###############################################################################
set -euo pipefail

REF_AMOA=${1:?"用法: bash run_oga_analysis_v3.sh /path/to/ref_all_amoA.fasta"}
THREADS=8
OUTDIR="./oga_results"
mkdir -p ${OUTDIR}

# ========== Step 1: BLASTn ==========
echo "Step 1: BLASTn分类..."
makeblastdb -in ${REF_AMOA} -dbtype nucl -out ${OUTDIR}/ref_db 2>/dev/null
blastn -query ABY89139.1_nucl_homologs.fna -db ${OUTDIR}/ref_db \
    -out ${OUTDIR}/blast_classification.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -max_target_seqs 1 -evalue 1e-10 -perc_identity 70 -num_threads ${THREADS}
echo "  BLASTn完成: $(wc -l < ${OUTDIR}/blast_classification.tsv) hits"

# ========== Step 2-4: Python ==========
python3 << 'PYEOF'
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import warnings, csv
warnings.filterwarnings('ignore')
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['font.size'] = 10

OUTDIR = "./oga_results"
COLORS = {'HighT':'#E8A07A','LowT':'#9CCFDE','NP_other':'#F4D03F','WCA':'#58D68D','WCB':'#8E44AD'}
ALL_ECOTYPES = ['HighT','LowT','NP_other','WCA','WCB']

# ===== 1. 读取BLASTn分类 =====
cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
        'qstart','qend','sstart','send','evalue','bitscore']
blast = pd.read_csv(f"{OUTDIR}/blast_classification.tsv", sep='\t', header=None, names=cols)
blast['ecotype'] = blast['sseqid'].str.split('__').str[0]
blast = blast.sort_values('bitscore', ascending=False).drop_duplicates('qseqid', keep='first')
classification = dict(zip(blast['qseqid'], blast['ecotype']))

print(f"分类结果: {len(classification)} unigenes")
for eco, n in blast['ecotype'].value_counts().items():
    print(f"  {eco}: {n}")
blast[['qseqid','sseqid','pident','ecotype']].to_csv(
    f"{OUTDIR}/classification_detail.tsv", sep='\t', index=False)

# ===== 2. 读取丰度矩阵 =====
with open("abundance_matrix.csv") as f:
    reader = csv.reader(f, delimiter='\t')
    row1 = next(reader)
    row2 = next(reader)
    
    sample_labels = row1[2:]
    sample_barcodes = row2[2:]
    
    data = {}
    for row in reader:
        if len(row) < 3:
            continue
        gene_id = row[0]
        vals = []
        for v in row[2:]:
            try:
                vals.append(float(v))
            except:
                vals.append(0.0)
        while len(vals) < len(sample_labels):
            vals.append(0.0)
        data[gene_id] = vals[:len(sample_labels)]

abund = pd.DataFrame.from_dict(data, orient='index', columns=sample_labels)
print(f"\n丰度矩阵: {abund.shape[0]} genes × {abund.shape[1]} samples")

# ===== 3. 只保留表层样本 (SRF + DCM) =====
surface_cols = [c for c in abund.columns if '_SRF_' in c or '_DCM_' in c]
abund_srf = abund[surface_cols].copy()
print(f"表层样本 (SRF+DCM): {len(surface_cols)}")

# ===== 4. 分类并汇总 =====
abund_srf['ecotype'] = abund_srf.index.map(classification)
n_class = abund_srf['ecotype'].notna().sum()
n_unclass = abund_srf['ecotype'].isna().sum()
print(f"已分类: {n_class}, 未分类: {n_unclass}")
abund_srf['ecotype'] = abund_srf['ecotype'].fillna('unclassified')

eco_sum = abund_srf.groupby('ecotype')[surface_cols].sum().T
eco_sum.index.name = 'sample'

# *** 修正: known包含NP_other ***
known = [e for e in ALL_ECOTYPES if e in eco_sum.columns]
eco_sum['total'] = eco_sum[known].sum(axis=1)

for eco in known:
    eco_sum[f'{eco}_frac'] = eco_sum[eco] / eco_sum['total'].replace(0, np.nan)

# *** 修正: NP = HighT + LowT + NP_other ***
np_cols = [e for e in ['HighT','LowT','NP_other'] if e in eco_sum.columns]
eco_sum['NP'] = eco_sum[np_cols].sum(axis=1)
eco_sum['NP_frac'] = eco_sum['NP'] / (eco_sum['NP'] + eco_sum.get('WCA', 0)).replace(0, np.nan)
eco_sum['WCA_frac_vs_NP'] = 1 - eco_sum['NP_frac']

# HighT within NP
ht = eco_sum.get('HighT', pd.Series(0, index=eco_sum.index))
lt = eco_sum.get('LowT', pd.Series(0, index=eco_sum.index))
eco_sum['HighT_in_NP'] = ht / (ht + lt).replace(0, np.nan)

eco_sum['log2_NP_WCA'] = np.log2((eco_sum['NP'] + 1e-12) / (eco_sum.get('WCA', 1e-12) + 1e-12))

# ===== 5. 合并metadata =====
env = pd.read_csv("environmental_parameters.csv", sep='\t')

barcode_to_env = {}
for _, row in env.iterrows():
    bid = str(row.get('sample_ID', ''))
    barcode_to_env[bid] = row

label_to_barcode = dict(zip(sample_labels, sample_barcodes))

for sample in eco_sum.index:
    barcode = label_to_barcode.get(sample, '')
    if barcode in barcode_to_env:
        row = barcode_to_env[barcode]
        eco_sum.loc[sample, 'latitude'] = row.get('latitude', np.nan)
        eco_sum.loc[sample, 'longitude'] = row.get('longitude', np.nan)
        eco_sum.loc[sample, 'temperature'] = row.get('Temperature (°C)', np.nan)
        eco_sum.loc[sample, 'depth_layer'] = row.get('depth', '')
        eco_sum.loc[sample, 'region'] = row.get('Region', '')
        eco_sum.loc[sample, 'dist_coast'] = row.get('Distance_coast (km)', np.nan)

# *** 修正: 去除重复样本 ***
eco_sum = eco_sum.reset_index()
eco_sum = eco_sum.drop_duplicates(subset='sample')
eco_sum = eco_sum.set_index('sample')

# 保存
eco_sum.to_csv(f"{OUTDIR}/surface_ecotype_abundance.tsv", sep='\t')

# ===== 6. 统计 =====
valid = eco_sum[eco_sum['total'] > 0].copy()
print(f"\n{'='*60}")
print(f"表层样本统计 ({len(valid)} 有效, 去重后)")
print(f"{'='*60}")

np_wca_valid = valid[valid['NP_frac'].notna()]
print(f"\nNP vs WCA (NP = HighT + LowT + NP_other):")
print(f"  NP fraction 中位数:  {np_wca_valid['NP_frac'].median():.4f}")
print(f"  WCA fraction 中位数: {np_wca_valid['WCA_frac_vs_NP'].median():.4f}")
wca_dom = (np_wca_valid['WCA_frac_vs_NP'] > 0.5).sum()
np_dom = (np_wca_valid['NP_frac'] > 0.5).sum()
print(f"  WCA占优: {wca_dom}/{len(np_wca_valid)} ({wca_dom/len(np_wca_valid)*100:.1f}%)")
print(f"  NP占优:  {np_dom}/{len(np_wca_valid)} ({np_dom/len(np_wca_valid)*100:.1f}%)")

print(f"\n各ecotype占总AOA reads中位数:")
for eco in known:
    fc = f'{eco}_frac'
    if fc in valid.columns:
        print(f"  {eco}: {valid[fc].median():.4f}")

# ===== 7. 出图 =====
print(f"\n生成图表...")

geo = valid.dropna(subset=['latitude','longitude']).copy()
geo['latitude'] = geo['latitude'].astype(float)
geo['longitude'] = geo['longitude'].astype(float)

# --- Fig 1: 全球地图（统一点大小）---
if len(geo) > 0:
    fig, ax = plt.subplots(figsize=(16, 8))
    cmap = LinearSegmentedColormap.from_list('c', ['#58D68D','#F5F5F5','#E8A07A'])
    
    sc = ax.scatter(geo['longitude'], geo['latitude'],
                   c=geo['NP_frac'], cmap=cmap, vmin=0, vmax=1,
                   s=70, edgecolors='#2E4053', linewidth=0.5, alpha=0.85)
    
    cbar = plt.colorbar(sc, ax=ax, shrink=0.6)
    cbar.set_label('NP fraction  (NP/(NP+WCA))')
    cbar.set_ticks([0,0.25,0.5,0.75,1])
    cbar.set_ticklabels(['WCA\ndominant','','Equal','','NP\ndominant'])
    
    ax.set_xlim(-180,180); ax.set_ylim(-75,80)
    ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
    ax.set_title(f'Nitrosopumilus vs WCA in global surface ocean (TARA, n={len(geo)})',
                fontsize=13, fontweight='bold')
    ax.set_facecolor('#f5f5f5')
    ax.axhline(0, color='gray', lw=0.3, ls=':')
    ax.text(0.02, 0.02,
            f'WCA dominant: {(geo["NP_frac"]<0.5).sum()}\nNP dominant: {(geo["NP_frac"]>0.5).sum()}',
            transform=ax.transAxes, fontsize=9, va='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    plt.tight_layout()
    for ext in ['pdf','png']:
        plt.savefig(f'{OUTDIR}/fig1_global_NP_WCA_surface.{ext}', dpi=300, bbox_inches='tight')
    print("  Fig 1: 全球地图 ✓")
    plt.close()

# --- Fig 2: 温度响应 ---
temp_df = geo.dropna(subset=['temperature']).copy()
temp_df['temperature'] = temp_df['temperature'].astype(float)

if len(temp_df) > 5:
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.scatter(temp_df['temperature'], temp_df['NP_frac'], c='#E8A07A', s=35, alpha=0.6,
              label='NP (HighT+LowT+NP_other)', edgecolors='none')
    ax.scatter(temp_df['temperature'], temp_df['WCA_frac_vs_NP'], c='#58D68D', s=35, alpha=0.6,
              label='WCA', edgecolors='none')
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        for col, color in [('NP_frac','#E8A07A'),('WCA_frac_vs_NP','#58D68D')]:
            s = temp_df.sort_values('temperature')
            sm = lowess(s[col].fillna(0), s['temperature'], frac=0.3)
            ax.plot(sm[:,0], sm[:,1], color=color, lw=2.5, zorder=10)
    except: pass
    ax.axhline(0.5, color='gray', lw=0.8, ls='--')
    ax.set_xlabel('Temperature (°C)', fontsize=12)
    ax.set_ylabel('Fraction (NP/(NP+WCA))', fontsize=12)
    ax.set_title('NP vs WCA across SST', fontweight='bold')
    ax.legend(frameon=False)
    ax.set_ylim(-0.05, 1.05)
    plt.tight_layout()
    for ext in ['pdf','png']:
        plt.savefig(f'{OUTDIR}/fig2_temperature_surface.{ext}', dpi=300, bbox_inches='tight')
    print("  Fig 2: 温度响应 ✓")
    plt.close()

# --- Fig 3: 纬度梯度（包含NP_other）---
lat_df = geo.copy()
if len(lat_df) > 5:
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios':[2,1]})
    
    ax = axes[0]
    for eco in ALL_ECOTYPES:
        fc = f'{eco}_frac'
        if fc in lat_df.columns:
            ax.scatter(lat_df['latitude'], lat_df[fc], c=COLORS.get(eco,'#ccc'),
                      s=20, alpha=0.5, label=eco, edgecolors='none')
            try:
                from statsmodels.nonparametric.smoothers_lowess import lowess
                s = lat_df.sort_values('latitude')
                sm = lowess(s[fc].fillna(0), s['latitude'], frac=0.3)
                ax.plot(sm[:,0], sm[:,1], color=COLORS.get(eco,'#ccc'), lw=2.5, zorder=10)
            except: pass
    ax.set_ylabel('Fraction of total AOA')
    ax.set_title('AOA ecotype latitude gradient (surface)')
    ax.legend(frameon=False, fontsize=8, ncol=3)
    
    ax2 = axes[1]
    ax2.scatter(lat_df['latitude'], lat_df['log2_NP_WCA'].clip(-8,8),
               c='#2E4053', s=15, alpha=0.5, edgecolors='none')
    ax2.axhline(0, color='gray', lw=1, ls='--')
    ax2.fill_between([-70,85], 0, 8, color='#E8A07A', alpha=0.04)
    ax2.fill_between([-70,85], -8, 0, color='#58D68D', alpha=0.04)
    ax2.set_xlabel('Latitude'); ax2.set_ylabel('log2(NP/WCA)')
    ax2.set_ylim(-8, 8)
    plt.tight_layout()
    for ext in ['pdf','png']:
        plt.savefig(f'{OUTDIR}/fig3_latitude_surface.{ext}', dpi=300, bbox_inches='tight')
    print("  Fig 3: 纬度梯度 ✓")
    plt.close()

# --- Fig 4: SRF vs DCM（包含NP_other）---
depth_df = valid.copy()
depth_df['layer'] = depth_df.index.map(lambda x: 'SRF' if '_SRF_' in x else 'DCM')

fig, ax = plt.subplots(figsize=(8, 5))
plot_ecos = [e for e in ALL_ECOTYPES if f'{e}_frac' in depth_df.columns]
x = np.arange(len(plot_ecos))
w = 0.35

for i, (layer, offset) in enumerate([('SRF', -w/2), ('DCM', w/2)]):
    sub = depth_df[depth_df['layer'] == layer]
    means = [sub[f'{eco}_frac'].median() if f'{eco}_frac' in sub.columns else 0 for eco in plot_ecos]
    colors = [COLORS.get(eco, '#ccc') for eco in plot_ecos]
    alpha = 0.85 if layer == 'SRF' else 0.55
    ax.bar(x + offset, means, w, label=layer, color=colors, alpha=alpha, edgecolor='white')

ax.set_xticks(x); ax.set_xticklabels(plot_ecos)
ax.set_ylabel('Median fraction'); ax.set_title('SRF vs DCM ecotype composition')
ax.legend(frameon=False)
plt.tight_layout()
for ext in ['pdf','png']:
    plt.savefig(f'{OUTDIR}/fig4_SRF_vs_DCM.{ext}', dpi=300, bbox_inches='tight')
print("  Fig 4: SRF vs DCM ✓")
plt.close()

print(f"\n所有结果在: {OUTDIR}/")
print("完成!")
PYEOF

echo ""
echo "=========================================="
echo "完成！结果在 ${OUTDIR}/"
echo "=========================================="
