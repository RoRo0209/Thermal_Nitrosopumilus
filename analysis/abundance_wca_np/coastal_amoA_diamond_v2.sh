#!/bin/bash
: "${THERMAL_PROJECT_ROOT:?Set THERMAL_PROJECT_ROOT to the analysis workspace}"
cd "${THERMAL_PROJECT_ROOT}/abundance_wca_np"

###############################################################################
# 近海数据 amoA DIAMOND BLASTx 竞争比对 v2
#
# 修正 (相对v1):
#   1. 参考库: ref_all_amoA.fasta (5组: HighT, LowT, NP_other, WCA, WCB)
#   2. 汇总包含 NP_other
#   3. NP = HighT + LowT + NP_other (与OGA v3一致)
#   4. 内置 total reads 计数 + percent归一化 (不再需要 add_normalization.sh)
#   5. Step 0 自动验证参考库分组完整性
#
# 用法:
#   bash coastal_amoA_diamond_v2.sh \
#     /path/to/ref_all_amoA.fasta \
#     /path/to/coastal_reads_dir \
#     [/path/to/sample_metadata.tsv]   # 可选，含 sample/latitude/longitude/temperature
#
# 参考库序列ID格式: {ecotype}__{seqname}
#   ecotype: HighT, LowT, NP_other, WCA, WCB
#
# 近海reads: *_1.fastq.gz / *_2.fastq.gz
###############################################################################

set -euo pipefail

REF_AMOA=${1:?"用法: bash coastal_amoA_diamond_v2.sh <ref_all_amoA.fasta> <reads_dir> [metadata.tsv]"}
READS_DIR=${2:?"请提供reads目录路径"}
METADATA=${3:-""}
THREADS=20
OUTDIR="./coastal_amoA_results_v2"

EVALUE="1e-10"
IDENTITY=80
R1_SUFFIX="_1.fastq.gz"
R2_SUFFIX="_2.fastq.gz"

mkdir -p ${OUTDIR}/{diamond_out,logs}

# ========== Step 0: 验证参考库 ==========
echo "============================================================"
echo "Step 0: 验证参考库分组"
echo "============================================================"
echo "  参考库: ${REF_AMOA}"
echo "  总序列数: $(grep -c '^>' ${REF_AMOA})"
echo "  各组:"

ALL_ECOTYPES="HighT LowT NP_other WCA WCB"
all_ok=true
for eco in ${ALL_ECOTYPES}; do
    n=$(grep -c "^>${eco}__" ${REF_AMOA} || true)
    echo "    ${eco}: ${n}"
    if [[ ${n} -eq 0 ]]; then
        echo "    *** WARNING: ${eco} 缺失! ***"
        all_ok=false
    fi
done
if [[ "${all_ok}" == "false" ]]; then
    echo "  参考库不完整，请检查!"
    exit 1
fi
echo "  OK"

# ========== Step 1: 翻译 + 建DIAMOND库 ==========
echo ""
echo "============================================================"
echo "Step 1: 准备amoA蛋白参考库"
echo "============================================================"

if [[ -f "${OUTDIR}/ref_amoA_db.dmnd" ]]; then
    echo "  DIAMOND数据库已存在，跳过"
else
    python3 << PYEOF
codon_table = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
def translate(seq):
    seq = seq.upper().replace('\n','')
    protein = []
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        aa = codon_table.get(codon, 'X')
        if aa == '*': break
        protein.append(aa)
    return ''.join(protein)

seqs, cur_id, cur_seq = {}, None, []
with open("${REF_AMOA}") as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if cur_id: seqs[cur_id] = ''.join(cur_seq)
            cur_id = line[1:].split()[0]; cur_seq = []
        else: cur_seq.append(line)
    if cur_id: seqs[cur_id] = ''.join(cur_seq)

n_out = 0
with open("${OUTDIR}/ref_amoA_protein.faa", 'w') as out:
    for sid, seq in seqs.items():
        prot = translate(seq)
        if len(prot) > 50:
            out.write(f">{sid}\n{prot}\n"); n_out += 1
print(f"翻译完成: {len(seqs)} 核酸 -> {n_out} 蛋白")
PYEOF

    diamond makedb \
        --in ${OUTDIR}/ref_amoA_protein.faa \
        --db ${OUTDIR}/ref_amoA_db \
        --threads ${THREADS} \
        2>${OUTDIR}/logs/makedb.log
    echo "  DIAMOND数据库构建完成"
fi

# ========== Step 2: DIAMOND BLASTx + 统计总reads ==========
echo ""
echo "============================================================"
echo "Step 2: DIAMOND BLASTx + 统计总reads"
echo "============================================================"

TOTAL_READS_FILE="${OUTDIR}/total_reads.tsv"
if [[ ! -f "${TOTAL_READS_FILE}" ]]; then
    echo -e "sample\ttotal_reads" > ${TOTAL_READS_FILE}
fi

sample_count=0
for r1 in ${READS_DIR}/*${R1_SUFFIX}; do
    [[ ! -f "$r1" ]] && continue
    base=$(basename "${r1}" "${R1_SUFFIX}")
    r2="${READS_DIR}/${base}${R2_SUFFIX}"

    # --- DIAMOND ---
    if [[ -f "${OUTDIR}/diamond_out/${base}.tsv" ]]; then
        echo "  跳过已完成: ${base}"
    else
        echo "  比对: ${base}"
        diamond blastx \
            --db ${OUTDIR}/ref_amoA_db \
            --query "${r1}" \
            --out ${OUTDIR}/diamond_out/${base}_R1.tsv \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --max-target-seqs 1 --evalue ${EVALUE} --id ${IDENTITY} \
            --threads ${THREADS} --sensitive \
            2>>${OUTDIR}/logs/diamond_${base}.log

        if [[ -f "${r2}" ]]; then
            diamond blastx \
                --db ${OUTDIR}/ref_amoA_db \
                --query "${r2}" \
                --out ${OUTDIR}/diamond_out/${base}_R2.tsv \
                --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
                --max-target-seqs 1 --evalue ${EVALUE} --id ${IDENTITY} \
                --threads ${THREADS} --sensitive \
                2>>${OUTDIR}/logs/diamond_${base}.log

            cat ${OUTDIR}/diamond_out/${base}_R1.tsv ${OUTDIR}/diamond_out/${base}_R2.tsv \
                > ${OUTDIR}/diamond_out/${base}.tsv
            rm -f ${OUTDIR}/diamond_out/${base}_R1.tsv ${OUTDIR}/diamond_out/${base}_R2.tsv
        else
            mv ${OUTDIR}/diamond_out/${base}_R1.tsv ${OUTDIR}/diamond_out/${base}.tsv
        fi
        hits=$(wc -l < ${OUTDIR}/diamond_out/${base}.tsv)
        echo "    hits: ${hits}"
    fi

    # --- 总reads数 ---
    if ! grep -q "^${base}	" ${TOTAL_READS_FILE} 2>/dev/null; then
        r1_count=$(zcat "${r1}" 2>/dev/null | awk 'END{print NR/4}')
        r2_count=0
        [[ -f "${r2}" ]] && r2_count=$(zcat "${r2}" 2>/dev/null | awk 'END{print NR/4}')
        total=$((r1_count + r2_count))
        echo -e "${base}\t${total}" >> ${TOTAL_READS_FILE}
    fi

    ((sample_count++)) || true
done
echo "  共处理 ${sample_count} 个样本"

# ========== Step 3: 汇总 (5组) + 归一化 ==========
echo ""
echo "============================================================"
echo "Step 3: 汇总 (5组含NP_other) + 归一化"
echo "============================================================"

METADATA_ARG="${METADATA}"

python3 << PYEOF
import pandas as pd
import numpy as np
import glob, os, sys

OUTDIR = "${OUTDIR}"
METADATA = "${METADATA_ARG}"
ALL_ECOTYPES = ['HighT', 'LowT', 'NP_other', 'WCA', 'WCB']

files = sorted(glob.glob(f"{OUTDIR}/diamond_out/*.tsv"))
files = [f for f in files if not f.endswith("_R1.tsv") and not f.endswith("_R2.tsv")]

cols = ['qseqid','sseqid','pident','length','mismatch','gapopen',
        'qstart','qend','sstart','send','evalue','bitscore']

results = []
for f in files:
    sample = os.path.basename(f).replace(".tsv", "")
    try:
        df = pd.read_csv(f, sep='\t', header=None, names=cols)
    except:
        continue
    if len(df) == 0:
        results.append({'sample': sample, 'total_amoA_reads': 0})
        continue

    # 提取ecotype: ID格式 ecotype__seqname
    # NP_other__xxx split by '__' -> ['NP_other', 'xxx'] -> [0] = 'NP_other' ✓
    df['ecotype'] = df['sseqid'].apply(lambda x: x.split('__')[0])

    # 检查未知ecotype
    unknown = set(df['ecotype'].unique()) - set(ALL_ECOTYPES)
    if unknown:
        print(f"  WARNING [{sample}]: 未知ecotype标签: {unknown}", file=sys.stderr)

    eco_counts = df['ecotype'].value_counts().to_dict()
    row = {'sample': sample, 'total_amoA_reads': len(df)}

    # *** 统计全部5组 ***
    for eco in ALL_ECOTYPES:
        row[eco] = eco_counts.get(eco, 0)

    # *** NP = HighT + LowT + NP_other (与OGA v3一致) ***
    row['NP_total'] = row['HighT'] + row['LowT'] + row['NP_other']
    row['WCA_total'] = row['WCA']
    row['WCB_total'] = row['WCB']

    # NP fraction (NP vs WCA, 不含WCB)
    np_wca = row['NP_total'] + row['WCA_total']
    row['NP_frac'] = row['NP_total'] / np_wca if np_wca > 0 else np.nan
    row['WCA_frac'] = row['WCA_total'] / np_wca if np_wca > 0 else np.nan

    # HighT within NP (只看HighT vs LowT, NP_other不参与)
    ht_lt = row['HighT'] + row['LowT']
    row['HighT_in_NP'] = row['HighT'] / ht_lt if ht_lt > 0 else np.nan

    # 各ecotype占total AOA
    total_known = sum(row[e] for e in ALL_ECOTYPES)
    for eco in ALL_ECOTYPES:
        row[f'{eco}_frac'] = row[eco] / total_known if total_known > 0 else np.nan

    row['data_source'] = 'coastal_DIAMOND'
    results.append(row)

coastal = pd.DataFrame(results)

# --- 合并总reads数 ---
try:
    tr = pd.read_csv(f"{OUTDIR}/total_reads.tsv", sep='\t')
    coastal = coastal.merge(tr[['sample','total_reads']], on='sample', how='left')
except:
    coastal['total_reads'] = np.nan

# --- 归一化 ---
if 'total_reads' in coastal.columns:
    coastal['total_amoA_abundance'] = coastal['total_amoA_reads'] / coastal['total_reads']
    for eco in ALL_ECOTYPES:
        coastal[f'{eco}_abundance'] = coastal[eco] / coastal['total_reads']
    coastal['NP_abundance'] = coastal['NP_total'] / coastal['total_reads']
    coastal['WCA_abundance'] = coastal['WCA_total'] / coastal['total_reads']

# --- metadata ---
if METADATA and os.path.exists(METADATA):
    meta = pd.read_csv(METADATA, sep='\t')
    if 'sample' in meta.columns:
        coastal = coastal.merge(meta, on='sample', how='left')
    print(f"  已合并metadata: {METADATA}")

# --- 保存 ---
outfile = f"{OUTDIR}/coastal_ecotype_summary.tsv"
coastal.to_csv(outfile, sep='\t', index=False)

# --- 统计 ---
print(f"\n{'='*60}")
print(f"近海样本汇总 (v2, 含NP_other修正)")
print(f"{'='*60}")
print(f"总样本数: {len(coastal)}")

valid = coastal[coastal['total_amoA_reads'] > 5]
print(f"有效样本 (>5 amoA hits): {len(valid)}")

if len(valid) > 0:
    print(f"\n各ecotype reads数中位数:")
    for eco in ALL_ECOTYPES:
        print(f"  {eco}: {valid[eco].median():.0f}")
    print(f"  NP_total: {valid['NP_total'].median():.0f}")

    print(f"\nNP vs WCA (NP = HighT + LowT + NP_other):")
    v = valid[valid['NP_frac'].notna()]
    print(f"  NP fraction 中位数: {v['NP_frac'].median():.4f}")
    np_dom = (v['NP_frac'] > 0.5).sum()
    wca_dom = (v['NP_frac'] < 0.5).sum()
    print(f"  NP占优: {np_dom}/{len(v)} ({np_dom/len(v)*100:.1f}%)")
    print(f"  WCA占优: {wca_dom}/{len(v)} ({wca_dom/len(v)*100:.1f}%)")

    ht_v = valid[valid['HighT_in_NP'].notna()]
    if len(ht_v) > 0:
        print(f"\nHighT within NP:")
        print(f"  HighT fraction 中位数: {ht_v['HighT_in_NP'].median():.4f}")

    if valid['total_reads'].notna().any():
        vn = valid[valid['total_reads'].notna()]
        print(f"\n归一化丰度 (reads/total_reads):")
        print(f"  Total amoA: {vn['total_amoA_abundance'].median():.2e}")
        print(f"  NP: {vn['NP_abundance'].median():.2e}")

print(f"\n结果: {outfile}")
PYEOF

echo ""
echo "=========================================="
echo "完成！"
echo "  结果: ${OUTDIR}/coastal_ecotype_summary.tsv"
echo "  NP = HighT + LowT + NP_other (与OGA v3一致)"
echo "  已含归一化 (不需要 add_normalization.sh)"
echo "=========================================="
