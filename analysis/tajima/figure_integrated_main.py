#!/usr/bin/env python
"""
figure_integrated_main.py
====================================================================
ISME J 整合主图 — 参考 SCBG 文献图1风格重绘

Layout:
  Row 1:       A. Manhattan plot (全基因组 Fst, 标注 C1-C6)
  Row 2:       B. π violin  |  C. SNV distance violin
  Row 3-4:     D. 6 个 differentiation island (3 col × 2 row)
               每个: 上半 Fst bar track + 下半基因箭头图
====================================================================
"""

# Release path adapter: data and outputs live in a separate user workspace.
import os as _release_os
from pathlib import Path as _ReleasePath
if not _release_os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the analysis workspace")
_release_root = _ReleasePath(_release_os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
_release_os.chdir(_release_root / 'tajima')


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from pathlib import Path
from itertools import combinations

Path("figures").mkdir(exist_ok=True)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 9,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
})

# ══════════════════════════════════════════════════════════════
# 配色
# ══════════════════════════════════════════════════════════════
PAL = {'LowT': '#5DADE2', 'HighT': '#E74C3C'}

FUNC_COLORS = {
    'Ammonia oxidation':     '#E74C3C',
    'Translation':           '#3498DB',
    'Transcription':         '#1F618D',
    'DNA replication/repair': '#27AE60',
    'Amino acid synthesis':   '#F39C12',
    'Cofactor synthesis':     '#8E44AD',
    'Transport':              '#1ABC9C',
    'Stress/folding':         '#E67E22',
    'Metabolism (other)':     '#7F8C8D',
    'Other/unknown':          '#B0B0B0',
    'Hypothetical':           '#D0D0D0',
}

CLUSTER_ACCENT = ['#E74C3C', '#E67E22', '#F1C40F', '#27AE60', '#3498DB', '#8E44AD']

# ══════════════════════════════════════════════════════════════
# 功能分类
# ══════════════════════════════════════════════════════════════
def classify(desc):
    if pd.isna(desc):
        return 'Hypothetical'
    d = desc.lower()
    if any(k in d for k in ['monooxygenase', 'amoa', 'amob', 'amoc', 'ammoni']):
        return 'Ammonia oxidation'
    if any(k in d for k in ['ribos', 'trna', 'rrna', 'translat', 'aminoacyl',
                            'elongat factor', 'initiation factor', 'elongation factor']):
        return 'Translation'
    if any(k in d for k in ['rna polymerase', 'nusa', 'spt5', 'transcription elongat']):
        return 'Transcription'
    if any(k in d for k in ['dna ligase', 'dna repair', 'replicat', 'helicase',
                            'methyltransferase', 'methylase', 'topoisomerase',
                            'rada', 'recombination', 'dna adenine']):
        return 'DNA replication/repair'
    if any(k in d for k in ['lysw', 'aminoadip', 'lysine', 'ornithine',
                            'argininosuccinate', 'hypusine', 'deoxyhypusine',
                            'carrier protein lysw']):
        return 'Amino acid synthesis'
    if any(k in d for k in ['riboflavin', 'lumazine', 'pyrimidinone',
                            'ribonucleoside-diphosphate reductase',
                            'cobal', 'b12', 'adenosylmethionine', 'sam synth',
                            'cyclohydrolase', 'butanone']):
        return 'Cofactor synthesis'
    if any(k in d for k in ['transport', 'permease', 'abc-2', 'mntr', 'manganese']):
        return 'Transport'
    if any(k in d for k in ['superoxide', 'prefoldin', 'chaperon']):
        return 'Stress/folding'
    if any(k in d for k in ['dehydrogenase', 'reductase', 'kinase', 'synthase',
                            'transferase', 'hydrolase', 'aldolase', 'lyase',
                            'carboxylase', 'carboxykinase', 'oxidase', 'desulfurase']):
        return 'Metabolism (other)'
    return 'Other/unknown'

# ══════════════════════════════════════════════════════════════
# 读数据
# ══════════════════════════════════════════════════════════════
fst = pd.read_csv("results/fst_windowed.tsv", sep="\t").dropna(subset=['Fst']).copy()
fst['Fst'] = fst['Fst'].clip(lower=0)
fst['mid'] = (fst['start'] + fst['end']) / 2

windowed = pd.read_csv("results/popgen_windowed.tsv", sep="\t")
outlier_genes = pd.read_csv("results/fst_outlier_genes_top5.tsv", sep="\t")
dist_df = pd.read_csv("results/pairwise_snv_distance.tsv", sep="\t", index_col=0)

manifest = pd.read_csv("sample_manifest.tsv", sep="\t")
highT_samples = set(manifest[manifest['clade'] == 'HighT']['sample_id']) & set(dist_df.index)
lowT_samples = set(manifest[manifest['clade'] == 'LowT']['sample_id']) & set(dist_df.index)

# GFF
gff_rows = []
with open('refs/LowT_representative.gff') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip().split('\t')
        if len(parts) < 9 or parts[2] != 'CDS':
            continue
        attrs = {}
        for kv in parts[8].split(';'):
            if '=' in kv:
                k, v = kv.split('=', 1)
                attrs[k.strip()] = v.strip()
        gff_rows.append({
            'gene_id': attrs.get('ID', ''),
            'strand': parts[6],
            'gene_name': attrs.get('Name', ''),
            'gff_start': int(parts[3]),
            'gff_end': int(parts[4]),
        })
gff_df = pd.DataFrame(gff_rows)

outlier_genes = outlier_genes.merge(
    gff_df[['gene_id', 'strand', 'gene_name']], on='gene_id', how='left')
outlier_genes['func_cat'] = outlier_genes['annotation'].apply(classify)

thr = fst['Fst'].quantile(0.95)

# ══════════════════════════════════════════════════════════════
# Cluster 定义
# ══════════════════════════════════════════════════════════════
# Study-selected intervals are an explicit input, not embedded results.
_clusters = pd.read_csv(_release_os.environ["THERMAL_FST_CLUSTERS"], sep="\t")
_required_cluster = {"id", "start", "end", "label"}
if not _required_cluster.issubset(_clusters.columns) or len(_clusters) != 6:
    raise ValueError("THERMAL_FST_CLUSTERS must contain six rows with id,start,end,label")
if _clusters["id"].tolist() != list(range(1, 7)) or not (_clusters["end"] >= _clusters["start"]).all():
    raise ValueError("Cluster IDs must be 1..6 in order and coordinates must be ordered")
CLUSTERS = _clusters.to_dict("records")

def assign_cluster(row):
    for cl in CLUSTERS:
        if row['gene_start'] <= cl['end'] + 5000 and row['gene_end'] >= cl['start']:
            return cl['id']
    return 0
outlier_genes['cluster'] = outlier_genes.apply(assign_cluster, axis=1)

# 关键基因标签 (只标这些，避免拥挤)
_gene_labels = pd.read_csv(_release_os.environ["THERMAL_FST_GENE_LABELS"], sep="\t")
if not {"gene_id", "label"}.issubset(_gene_labels.columns) or _gene_labels["gene_id"].duplicated().any():
    raise ValueError("THERMAL_FST_GENE_LABELS requires unique gene_id and label columns")
KEY_GENES = dict(zip(_gene_labels["gene_id"], _gene_labels["label"]))

# ══════════════════════════════════════════════════════════════
# 画基因箭头
# ══════════════════════════════════════════════════════════════
def draw_gene_arrow(ax, x_start, x_end, y, height, strand, color,
                    edgecolor='#333333', linewidth=0.6, alpha=0.9):
    width = x_end - x_start
    if width <= 0:
        return
    arrow_head = min(width * 0.2, 1200)

    if strand == '+':
        body_end = x_start + width - arrow_head
        verts = [
            (x_start, y - height / 2),
            (max(body_end, x_start), y - height / 2),
            (x_end, y),
            (max(body_end, x_start), y + height / 2),
            (x_start, y + height / 2),
        ]
    else:
        body_start = x_start + arrow_head
        verts = [
            (x_end, y - height / 2),
            (min(body_start, x_end), y - height / 2),
            (x_start, y),
            (min(body_start, x_end), y + height / 2),
            (x_end, y + height / 2),
        ]

    poly = Polygon(verts, closed=True, facecolor=color, edgecolor=edgecolor,
                   linewidth=linewidth, alpha=alpha, zorder=3)
    ax.add_patch(poly)


# ══════════════════════════════════════════════════════════════
# 图布局
# ══════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(18, 36))

gs = GridSpec(8, 2, figure=fig,
              height_ratios=[1.2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
              hspace=0.35, wspace=0.25,
              left=0.06, right=0.94, top=0.97, bottom=0.03)

ax_A = fig.add_subplot(gs[0, :])       # Manhattan 占满第一行
ax_B = fig.add_subplot(gs[1, 0])       # π
ax_C = fig.add_subplot(gs[1, 1])       # SNV distance
# 6 cluster subplots (1 col × 6 rows, 每个占满整行宽度)
ax_cl = []
for i in range(6):
    ax_cl.append(fig.add_subplot(gs[2 + i, :]))

# ══════════════════════════════════════════════════════════════
# Panel A: Manhattan
# ══════════════════════════════════════════════════════════════
ax = ax_A
fst_s = fst.sort_values('mid').copy()

# 背景点
bg = fst_s[fst_s['Fst'] < thr]
ax.scatter(bg['mid'], bg['Fst'], color='#BDC3C7', s=18, alpha=0.6,
           edgecolor='none', zorder=1)

# outlier 点
ol = fst_s[fst_s['Fst'] >= thr]
ax.scatter(ol['mid'], ol['Fst'], color='#E74C3C', s=30, alpha=0.9,
           edgecolor='white', linewidth=0.4, zorder=2)

# 阈值线
ax.axhline(thr, color='#E74C3C', linestyle='--', linewidth=1, alpha=0.6)
ax.text(fst_s['mid'].max() * 0.99, thr + 0.008,
        f'Top 5% = {thr:.3f}', fontsize=8, ha='right', color='#E74C3C')

# 中位线
ax.axhline(fst_s['Fst'].median(), color='#7F8C8D', linestyle=':', linewidth=0.8, alpha=0.5)

# 标注 cluster
for i, cl in enumerate(CLUSTERS):
    ax.axvspan(cl['start'], cl['end'], alpha=0.12, color=CLUSTER_ACCENT[i], zorder=0)
    mid = (cl['start'] + cl['end']) / 2
    ax.annotate(f'C{cl["id"]}', xy=(mid, fst_s['Fst'].max() + 0.01),
                fontsize=9, fontweight='bold', ha='center', va='bottom',
                color=CLUSTER_ACCENT[i])

ax.set_xlim(0, fst_s['mid'].max() * 1.02)
ax.set_ylim(-0.02, fst_s['Fst'].max() + 0.05)
xticks = np.arange(0, fst_s['mid'].max(), 200000)
ax.set_xticks(xticks)
ax.set_xticklabels([f'{int(x/1000)}' for x in xticks])
ax.set_xlabel('Genome position (kb)')
ax.set_ylabel('Fst (HighT vs LowT)')
ax.set_title('A', fontweight='bold', fontsize=13, loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ══════════════════════════════════════════════════════════════
# Panel B: π violin
# ══════════════════════════════════════════════════════════════
ax = ax_B
for idx, (clade, pos) in enumerate(zip(['LowT', 'HighT'], [0, 1])):
    d = windowed[(windowed.clade == clade) & windowed.pi.notna()]['pi']
    if len(d) == 0:
        continue
    parts = ax.violinplot(d, positions=[pos], widths=0.6,
                          showmedians=False, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor(PAL[clade])
        pc.set_alpha(0.7)
        pc.set_edgecolor('white')
        pc.set_linewidth(0.8)
    med = d.median()
    ax.hlines(med, pos - 0.18, pos + 0.18, color='black', linewidth=2.5, zorder=5)
    ax.text(pos, med + 0.0015, f'{med:.4f}', ha='center', fontsize=8,
            fontweight='bold', zorder=6)

ax.set_xticks([0, 1])
ax.set_xticklabels(['LowT\n(n=19)', 'HighT\n(n=7)'])
ax.set_ylabel('π (per 10-kb window)')
ax.set_title('B', fontweight='bold', fontsize=13, loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.text(0.97, 0.95,
        'Subsampling validated\n(100 replicates, p < 0.01)',
        transform=ax.transAxes, fontsize=7, va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FAFAFA',
                  edgecolor='#CCCCCC', alpha=0.9))

# ══════════════════════════════════════════════════════════════
# Panel C: SNV distance
# ══════════════════════════════════════════════════════════════
ax = ax_C
within_ht, within_lt, between_dist = [], [], []
ht_list, lt_list = sorted(highT_samples), sorted(lowT_samples)

for s1, s2 in combinations(ht_list, 2):
    if s1 in dist_df.index and s2 in dist_df.columns:
        within_ht.append(dist_df.loc[s1, s2])
for s1, s2 in combinations(lt_list, 2):
    if s1 in dist_df.index and s2 in dist_df.columns:
        within_lt.append(dist_df.loc[s1, s2])
for s1 in ht_list:
    for s2 in lt_list:
        if s1 in dist_df.index and s2 in dist_df.columns:
            between_dist.append(dist_df.loc[s1, s2])

groups = [within_ht, within_lt, between_dist]
labels_c = ['Within\nHighT', 'Within\nLowT', 'Between']
colors_c = [PAL['HighT'], PAL['LowT'], '#95A5A6']

for i, (data, lbl, col) in enumerate(zip(groups, labels_c, colors_c)):
    if not data:
        continue
    parts = ax.violinplot(data, positions=[i], widths=0.65,
                          showmedians=False, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor(col)
        pc.set_alpha(0.7)
        pc.set_edgecolor('white')
        pc.set_linewidth(0.8)
    med = np.median(data)
    ax.hlines(med, i - 0.18, i + 0.18, color='black', linewidth=2.5, zorder=5)
    ax.text(i, med + 0.006, f'{med:.1%}', ha='center', fontsize=8,
            fontweight='bold', zorder=6)

ax.set_xticks([0, 1, 2])
ax.set_xticklabels(labels_c)
ax.set_ylabel('Pairwise SNV distance')
ax.set_title('C', fontweight='bold', fontsize=13, loc='left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 显著性标注
def add_bracket(ax, x1, x2, y, h, text):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.8, color='black')
    ax.text((x1+x2)/2, y+h+0.002, text, ha='center', fontsize=7)

ymax = max(max(within_ht), max(within_lt), max(between_dist)) + 0.01
add_bracket(ax, 0, 2, ymax, 0.01, '***')
add_bracket(ax, 1, 2, ymax + 0.03, 0.01, 'n.s.')

# ══════════════════════════════════════════════════════════════
# Panel D: 6 个 cluster 基因簇图
# ══════════════════════════════════════════════════════════════
outlier_id_set = set(outlier_genes['gene_id'])

for ci, (cl, ax) in enumerate(zip(CLUSTERS, ax_cl)):
    pad = 3000
    region_start = cl['start'] - pad
    region_end = cl['end'] + pad

    # ── 上半: Fst bar track ──
    fst_region = fst_s[(fst_s['mid'] >= region_start - 5000) &
                       (fst_s['mid'] <= region_end + 5000)].copy()

    # 用 twinx 画 Fst bars
    ax_fst = ax.twinx()

    if len(fst_region) > 0:
        bar_w = (fst_region['mid'].diff().median() or 5000) * 0.8
        colors_bar = []
        for _, r in fst_region.iterrows():
            if r['Fst'] >= thr:
                colors_bar.append(CLUSTER_ACCENT[ci])
            else:
                colors_bar.append('#D5D8DC')
        ax_fst.bar(fst_region['mid'], fst_region['Fst'],
                   width=bar_w, color=colors_bar, edgecolor='white',
                   linewidth=0.3, alpha=0.7, zorder=1)
        ax_fst.axhline(thr, color='#E74C3C', linestyle='--',
                       linewidth=0.7, alpha=0.4)

    ax_fst.set_ylim(0, 0.50)
    ax_fst.set_ylabel('Fst', fontsize=7, color='#7F8C8D', rotation=0, labelpad=12)
    ax_fst.tick_params(axis='y', labelsize=6, colors='#7F8C8D')
    ax_fst.spines['top'].set_visible(False)

    # ── 下半: 基因箭头 ──
    all_genes = gff_df[(gff_df['gff_start'] <= region_end + 2000) &
                       (gff_df['gff_end'] >= region_start - 2000)].copy()

    # merge 功能注释
    all_genes = all_genes.merge(
        outlier_genes[['gene_id', 'func_cat', 'KO', 'annotation']].drop_duplicates('gene_id'),
        on='gene_id', how='left')
    all_genes['is_outlier'] = all_genes['gene_id'].isin(outlier_id_set)
    all_genes['func_cat'] = all_genes['func_cat'].fillna('Hypothetical')

    arrow_y = -0.35
    arrow_h = 0.28

    for _, gene in all_genes.sort_values('gff_start').iterrows():
        if gene['is_outlier']:
            color = FUNC_COLORS.get(gene['func_cat'], '#B0B0B0')
            ec = '#333333'
            lw = 0.7
        else:
            color = '#EEEEEE'
            ec = '#CCCCCC'
            lw = 0.4

        draw_gene_arrow(ax, gene['gff_start'], gene['gff_end'],
                        arrow_y, arrow_h, gene['strand'], color,
                        edgecolor=ec, linewidth=lw)

        # 标签: 只标关键基因
        label = KEY_GENES.get(gene['gene_id'], '')
        if label:
            mid_x = (gene['gff_start'] + gene['gff_end']) / 2
            ax.annotate(label,
                        xy=(mid_x, arrow_y + arrow_h / 2),
                        xytext=(0, 8), textcoords='offset points',
                        fontsize=6, ha='center', va='bottom',
                        fontstyle='italic', fontweight='bold',
                        color='#2C3E50',
                        arrowprops=dict(arrowstyle='-', color='#7F8C8D',
                                        linewidth=0.4, shrinkA=0, shrinkB=1))

    # 坐标轴
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-0.6, 0.15)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # x 轴 kb
    xt = np.linspace(region_start, region_end, 4)
    ax.set_xticks(xt)
    ax.set_xticklabels([f'{x/1000:.0f}' for x in xt], fontsize=7)
    ax.set_xlabel('Position (kb)', fontsize=7)

    # 标题
    peak_fst = fst_region['Fst'].max() if len(fst_region) > 0 else 0
    ax.set_title(f'C{cl["id"]}: {cl["label"]}',
                 fontweight='bold', fontsize=9, loc='left',
                 color=CLUSTER_ACCENT[ci], pad=8)

    # 基因组线
    ax.axhline(arrow_y - arrow_h / 2 - 0.02, color='black',
               linewidth=0.4, alpha=0.3, zorder=0)

# ══════════════════════════════════════════════════════════════
# 功能图例
# ══════════════════════════════════════════════════════════════
legend_cats = ['Ammonia oxidation', 'Translation', 'Transcription',
               'DNA replication/repair', 'Amino acid synthesis', 'Cofactor synthesis',
               'Transport', 'Stress/folding', 'Metabolism (other)', 'Other/unknown']
legend_elements = [mpatches.Patch(facecolor=FUNC_COLORS[k], edgecolor='#333333',
                                   linewidth=0.5, label=k)
                   for k in legend_cats]
legend_elements.append(mpatches.Patch(facecolor='#EEEEEE', edgecolor='#CCCCCC',
                                       linewidth=0.5, label='Non-outlier'))

fig.legend(handles=legend_elements, loc='lower center', ncol=6, fontsize=7.5,
           frameon=True, framealpha=0.95, edgecolor='#CCCCCC',
           bbox_to_anchor=(0.5, 0.005), handlelength=1.5, handleheight=1.0,
           columnspacing=1.0)

# ── 保存 ──
try:
    plt.savefig('figures/figure_integrated_main.pdf', bbox_inches='tight', dpi=300)
except PermissionError:
    print("  PDF locked, skipping")
plt.savefig('figures/figure_integrated_main.png', bbox_inches='tight', dpi=200)
plt.close()
print("Done → figures/figure_integrated_main.{pdf,png}")
