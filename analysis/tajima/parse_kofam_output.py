#!/usr/bin/env python
"""
parse_kofam_output.py
====================================================================
解析 KofamKOALA 输出，输出 clean KO annotation table

输入:
  - detail-tsv 格式: *.txt (带 * 标记高置信度)
  或
  - mapper-one-line 格式: *.txt (每行一个 gene + tab + KO)

输出:
  - {prefix}_kofam_clean.tsv  干净的 KO 表

用法:
  python parse_kofam_output.py --detail LowT_detail.txt --output LowT_kofam
  python parse_kofam_output.py --mapper LowT_mapper.txt --output LowT_kofam
  python parse_kofam_output.py --detail X.txt --mapper Y.txt --output LowT_kofam
====================================================================
"""

# Release path adapter: data and outputs live in a separate user workspace.
import os as _release_os
from pathlib import Path as _ReleasePath
if not _release_os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the analysis workspace")
_release_root = _ReleasePath(_release_os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
_release_os.chdir(_release_root / 'tajima')


import argparse
import pandas as pd
import sys
import re

parser = argparse.ArgumentParser()
parser.add_argument('--detail', help='KofamKOALA detail-tsv output (带 * 标记)')
parser.add_argument('--mapper', help='KofamKOALA mapper-one-line output')
parser.add_argument('--output', required=True, help='输出文件前缀')
parser.add_argument('--strip_prefix', default='',
                    help='从 gene_id 去掉的前缀 (常用于匹配 GFF)，支持 regex')
parser.add_argument('--only_high_conf', action='store_true', default=True,
                    help='只保留高置信度 (带 * 的行)')
args = parser.parse_args()

if not args.detail and not args.mapper:
    print("错误：至少提供 --detail 或 --mapper 之一")
    sys.exit(1)

# -------------------- 解析 detail 文件 --------------------
detail_df = None
if args.detail:
    print(f"解析 detail 文件: {args.detail}")
    rows = []
    with open(args.detail) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            is_hc = line.startswith("*")
            text = line[1:].strip() if is_hc else line.strip()
            parts = text.split(None, 5)
            if len(parts) < 6:
                continue
            try:
                rows.append({
                    "gene_id": parts[0],
                    "KO": parts[1],
                    "thrshld": float(parts[2]),
                    "score": float(parts[3]),
                    "evalue": parts[4],
                    "description": parts[5],
                    "high_conf": is_hc,
                })
            except ValueError:
                continue
    detail_df = pd.DataFrame(rows)
    print(f"  总 hits: {len(detail_df)}, 高置信度: {detail_df['high_conf'].sum()}")

# -------------------- 解析 mapper 文件 --------------------
mapper_df = None
if args.mapper:
    print(f"解析 mapper 文件: {args.mapper}")
    rows = []
    with open(args.mapper) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 1:
                rows.append({"gene_id": parts[0], "KO": None})
            elif len(parts) == 2:
                rows.append({"gene_id": parts[0], "KO": parts[1]})
    mapper_df = pd.DataFrame(rows)
    print(f"  总 gene 数: {len(mapper_df)}, 有 KO: {mapper_df['KO'].notna().sum()}")

# -------------------- 合并策略 --------------------
# 如果两者都有: 用 mapper 作为全部 gene list (包括没 KO 的)，用 detail 补充 description 和 score
# 如果只有 detail: 直接用 detail
# 如果只有 mapper: 只有 KO 没 description

if detail_df is not None and mapper_df is not None:
    # 先只保留 detail 里高置信度
    if args.only_high_conf:
        detail_use = detail_df[detail_df['high_conf']]
    else:
        detail_use = detail_df
    # 按 gene_id 聚合（一个 gene 可能有多个 KO hit，取 score 最高的）
    detail_best = detail_use.sort_values('score', ascending=False).drop_duplicates('gene_id', keep='first')
    # merge to mapper
    merged = mapper_df.merge(
        detail_best[['gene_id', 'KO', 'description', 'score', 'evalue']],
        on='gene_id', how='left', suffixes=('_mapper', '_detail')
    )
    # 优先用 detail 的 KO
    merged['KO'] = merged['KO_detail'].fillna(merged['KO_mapper'])
    # 清理 KO_mapper/KO_detail
    final = merged[['gene_id', 'KO', 'description', 'score', 'evalue']].copy()
elif detail_df is not None:
    if args.only_high_conf:
        detail_use = detail_df[detail_df['high_conf']]
    else:
        detail_use = detail_df
    final = detail_use.sort_values('score', ascending=False).drop_duplicates('gene_id', keep='first')
    final = final[['gene_id', 'KO', 'description', 'score', 'evalue']].copy()
else:
    final = mapper_df.copy()
    final['description'] = None
    final['score'] = None
    final['evalue'] = None

# -------------------- 处理 gene_id 前缀 --------------------
# 提取 gene symbol (description 里首个词，比如 "amoA; ammonia...")
def extract_symbol(desc):
    if pd.isna(desc):
        return None
    # 匹配 xxx; 格式的 gene symbol
    m = re.match(r'^([a-zA-Z][\w\d\-]+);', desc)
    if m:
        return m.group(1)
    # 或者取第一个空格前的词（小心不要取到全大写的 "NAD" 之类）
    return None

final['symbol'] = final['description'].apply(extract_symbol)

# 如果需要去掉前缀（比如 kofam 输出是 "GCF_xxx__PREFIX_00001"，GFF 里只有 "PREFIX_00001"）
if args.strip_prefix:
    print(f"从 gene_id 去掉前缀: {args.strip_prefix}")
    final['gene_id_original'] = final['gene_id']
    final['gene_id'] = final['gene_id'].str.replace(args.strip_prefix, '', regex=True)

# -------------------- 输出 --------------------
out_file = f"{args.output}_kofam_clean.tsv"
final.to_csv(out_file, sep='\t', index=False)
print(f"\n输出: {out_file}")
print(f"  总 gene 数:        {len(final)}")
print(f"  有 KO 注释:        {final['KO'].notna().sum()}")
print(f"  注释率:            {final['KO'].notna().mean()*100:.1f}%")

# -------------------- 展示前 5 行 --------------------
print(f"\n前 5 行:")
print(final.head().to_string(index=False))

# -------------------- 额外：KO/pathway 简单统计 --------------------
print(f"\n最常见的 10 个 KO:")
top_ko = final[final['KO'].notna()].groupby(['KO', 'description']).size().sort_values(ascending=False).head(10)
for (ko, desc), n in top_ko.items():
    print(f"  {ko}  ×{n}  {desc[:70]}")
