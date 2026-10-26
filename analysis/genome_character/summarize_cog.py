#!/usr/bin/env python3
"""
summarize_cog.py

从 eggNOG-mapper 的注释结果中，统计每个基因组的 COG 类别基因数量，
输出用于 Y-A-S 生活史策略分析的汇总表。

用法:
    python3 summarize_cog.py eggnog_results/

输出:
    cog_summary_for_YAS.tsv   — 每个基因组 × 每个COG类别的基因计数
    yas_scores.tsv             — 每个基因组的S-strategy相关指标 (COG_T, COG_L, COG_O)
                                  以及其他所有COG类别，方便你合并到总表

COG 功能类别说明 (与 Y-A-S 的对应关系):
    Y-strategy 相关:
        J - Translation, ribosomal structure and biogenesis
        (配合 genome_size, coding_density, gRodon 等已有指标)
    A-strategy 相关:
        G - Carbohydrate transport and metabolism
        E - Amino acid transport and metabolism
        P - Inorganic ion transport and metabolism
        (配合 transporter, CAZyme, peptidase 等已有指标)
    S-strategy 相关:
        T - Signal transduction mechanisms
        L - Replication, recombination and repair
        O - Post-translational modification, protein turnover, chaperones
"""

import os
from pathlib import Path
_thermal_root = os.environ.get("THERMAL_PROJECT_ROOT", "")
if not _thermal_root or not Path(_thermal_root).is_dir():
    raise SystemExit("Set THERMAL_PROJECT_ROOT to an existing, separate analysis workspace.")
THERMAL_PROJECT_ROOT = Path(_thermal_root).resolve()

import sys
import csv
from collections import defaultdict

# 所有 COG 功能类别
os.chdir(THERMAL_PROJECT_ROOT / 'genome_character')

COG_CATEGORIES = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

COG_DESCRIPTIONS = {
    'A': 'RNA processing and modification',
    'B': 'Chromatin structure and dynamics',
    'C': 'Energy production and conversion',
    'D': 'Cell cycle control, cell division',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'G': 'Carbohydrate transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'J': 'Translation, ribosomal structure',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',  # S-strategy
    'M': 'Cell wall/membrane biogenesis',
    'N': 'Cell motility',
    'O': 'PTM, protein turnover, chaperones',      # S-strategy
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'T': 'Signal transduction mechanisms',          # S-strategy
    'U': 'Intracellular trafficking, secretion',
    'V': 'Defense mechanisms',
    'W': 'Extracellular structures',
    'X': 'Mobilome: prophages, transposons',
    'Y': 'Nuclear structure',
    'Z': 'Cytoskeleton',
}


def parse_eggnog_annotations(filepath):
    """
    解析 eggNOG-mapper 的 .emapper.annotations 文件
    返回每个 COG 类别的基因计数
    """
    cog_counts = defaultdict(int)
    total_genes = 0
    annotated_genes = 0
    
    with open(filepath, 'r') as f:
        for line in f:
            # 跳过注释行
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 7:
                continue
            
            total_genes += 1
            
            # eggNOG-mapper v2 输出格式:
            # 第7列 (index 6) 是 COG_category
            # 但不同版本列号可能不同，我们用列名来定位
            # 通常是: query | seed_ortholog | evalue | score | eggNOG_OGs | 
            #         max_annot_lvl | COG_category | Description | ...
            
            cog_cat = fields[6].strip() if len(fields) > 6 else '-'
            
            if cog_cat and cog_cat != '-' and cog_cat != '':
                annotated_genes += 1
                # 一个基因可能属于多个 COG 类别 (如 "GT" 表示同时属于 G 和 T)
                for c in cog_cat:
                    if c in COG_CATEGORIES:
                        cog_counts[c] += 1
    
    return dict(cog_counts), total_genes, annotated_genes


def main():
    if len(sys.argv) < 2:
        print("用法: python3 summarize_cog.py <eggnog_results_dir>")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    
    # 找到所有 .emapper.annotations 文件
    ann_files = []
    for f in sorted(os.listdir(results_dir)):
        if f.endswith('.emapper.annotations'):
            ann_files.append(f)
    
    if not ann_files:
        print(f"错误: 在 {results_dir} 中没找到 .emapper.annotations 文件")
        sys.exit(1)
    
    print(f"找到 {len(ann_files)} 个注释文件")
    
    # 解析每个文件
    all_data = {}
    for ann_file in ann_files:
        # 从文件名提取基因组名
        # 例如: SCM1.emapper.annotations -> SCM1
        genome_name = ann_file.replace('.emapper.annotations', '')
        
        filepath = os.path.join(results_dir, ann_file)
        cog_counts, total_genes, annotated_genes = parse_eggnog_annotations(filepath)
        
        all_data[genome_name] = {
            'cog_counts': cog_counts,
            'total_genes': total_genes,
            'annotated_genes': annotated_genes,
        }
        
        print(f"  {genome_name}: {total_genes} genes, {annotated_genes} annotated, "
              f"COG_T={cog_counts.get('T', 0)}, COG_L={cog_counts.get('L', 0)}, "
              f"COG_O={cog_counts.get('O', 0)}")
    
    # ============================================================
    # 输出1: 完整的 COG 类别计数表
    # ============================================================
    out1 = 'cog_summary_for_YAS.tsv'
    with open(out1, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        header = ['Genome', 'Total_genes', 'Annotated_genes'] + \
                 [f'COG_{c}' for c in COG_CATEGORIES]
        writer.writerow(header)
        
        for genome_name in sorted(all_data.keys()):
            data = all_data[genome_name]
            row = [genome_name, data['total_genes'], data['annotated_genes']]
            for c in COG_CATEGORIES:
                row.append(data['cog_counts'].get(c, 0))
            writer.writerow(row)
    
    print(f"\n完整COG统计表已保存: {out1}")
    
    # ============================================================
    # 输出2: S-strategy 相关指标提取
    # ============================================================
    out2 = 'yas_s_strategy_indicators.tsv'
    with open(out2, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'Genome',
            'Total_genes',
            'COG_T_signal_transduction',
            'COG_L_replication_repair',
            'COG_O_PTM_chaperones',
            'COG_V_defense',
            'COG_N_motility',
            'COG_K_transcription',
            'S_strategy_sum_TLO',
        ])
        
        for genome_name in sorted(all_data.keys()):
            cc = all_data[genome_name]['cog_counts']
            t_count = cc.get('T', 0)
            l_count = cc.get('L', 0)
            o_count = cc.get('O', 0)
            v_count = cc.get('V', 0)
            n_count = cc.get('N', 0)
            k_count = cc.get('K', 0)
            
            writer.writerow([
                genome_name,
                all_data[genome_name]['total_genes'],
                t_count,
                l_count,
                o_count,
                v_count,
                n_count,
                k_count,
                t_count + l_count + o_count,
            ])
    
    print(f"S-strategy指标表已保存: {out2}")
    
    # ============================================================
    # These counts are inputs to a separately specified trait analysis.


if __name__ == '__main__':
    main()
