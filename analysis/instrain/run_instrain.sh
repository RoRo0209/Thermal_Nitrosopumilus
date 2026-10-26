#!/bin/bash
: "${THERMAL_PROJECT_ROOT:?Set THERMAL_PROJECT_ROOT to the analysis workspace}"
cd "${THERMAL_PROJECT_ROOT}/instrain"

# ==============================================================================
# run_instrain.sh
# ==============================================================================
#
# 26 个样本（和 Fst 一致，来自 vcf_shared/HighT_LowT_joint.vcf.gz）
# ref:  "${THERMAL_AOA_REFERENCE_FASTA:?Set THERMAL_AOA_REFERENCE_FASTA to the required input}"（全部 HighT + LowT 基因组拼成）
# stb:  "${THERMAL_AOA_STB:?Set THERMAL_AOA_STB to the required input}"
# gene: "${THERMAL_AOA_GENE_FASTA:?Set THERMAL_AOA_GENE_FASTA to the required input}"
#
# 流程:
#   Step 1.  把 26 个样本的 reads 重新 map 到 all_AOA_contig.fa
#            （Fst 用的 BAM 是 map 到单个 LowT ref 的，不能直接用）
#   Step 2.  inStrain profile per-sample
# ==============================================================================

set -uo pipefail

# -------------------- 配置 --------------------
SAMPLES_TSV="fst_26samples.tsv"
MANIFEST="sample_manifest.tsv"     # 用来查 reads 路径

REF="${THERMAL_AOA_REFERENCE_FASTA:?Set THERMAL_AOA_REFERENCE_FASTA to the required input}"
STB="${THERMAL_AOA_STB:?Set THERMAL_AOA_STB to the required input}"
GENES="${THERMAL_AOA_GENE_FASTA:?Set THERMAL_AOA_GENE_FASTA to the required input}"

REF_BASE="${REF%.fa}"
THREADS=36
MIN_MAPQ=2          # inStrain 推荐 -q 2（保留 multi-mapped 但 quality 过滤）
                    # 也可以更严：-q 20

BAM_DIR="bam_instrain"          # map 到 all_AOA_contig 的 BAM
PROFILE_DIR="instrain_profile"
LOG_DIR="logs_instrain"

mkdir -p $BAM_DIR $PROFILE_DIR $LOG_DIR results

# -------------------- 检查依赖 --------------------
echo "[$(date +%H:%M:%S)] 检查依赖..."
for cmd in inStrain bowtie2 samtools; do
    command -v $cmd >/dev/null 2>&1 || { echo "  ⚠ 缺少 $cmd"; exit 1; }
done
for f in $REF $STB $GENES $SAMPLES_TSV $MANIFEST; do
    [ -f "$f" ] || { echo "  ⚠ 缺少文件: $f"; exit 1; }
done
echo "  OK"

# -------------------- Step 1: bowtie2 索引 --------------------
if [ ! -f ${REF_BASE}.1.bt2 ] && [ ! -f ${REF_BASE}.1.bt2l ]; then
    echo ""
    echo "[$(date +%H:%M:%S)] 建 bowtie2 索引 (这步对 all_AOA_contig 是大文件，可能要几十分钟)..."
    bowtie2-build --threads $THREADS $REF $REF_BASE
fi
[ -f ${REF}.fai ] || samtools faidx $REF

# -------------------- Step 2: per-sample mapping --------------------
echo ""
echo "[$(date +%H:%M:%S)] 重新 mapping 26 个样本到 all_AOA_contig.fa..."

# 用 awk 从 manifest 查 reads 路径
get_reads() {
    local s=$1
    awk -F'\t' -v s=$s 'NR>1 && $1==s {print $2"\t"$3; exit}' $MANIFEST
}

while IFS=$'\t' read -r sample clade rest; do
    [[ $sample == "sample_id" ]] && continue
    
    BAM=${BAM_DIR}/${sample}.bam
    if [ -f $BAM ] && [ -f ${BAM}.bai ]; then
        echo "  [skip] $sample"
        continue
    fi
    
    reads=$(get_reads $sample)
    r1=$(echo "$reads" | cut -f1)
    r2=$(echo "$reads" | cut -f2)
    
    if [ -z "$r1" ] || [ -z "$r2" ]; then
        echo "  ⚠ $sample: reads 路径找不到"
        continue
    fi
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        echo "  ⚠ $sample: reads 文件不存在"
        continue
    fi
    
    echo "  [map] $sample ($clade)"
    bowtie2 -x $REF_BASE \
        -1 $r1 -2 $r2 \
        --very-sensitive --no-unal \
        -p $THREADS 2> ${LOG_DIR}/${sample}_bowtie.log | \
    samtools view -b -q $MIN_MAPQ -F 4 -@ 4 - | \
    samtools sort -@ 4 -o ${BAM}.tmp -
    
    pipe_status=("${PIPESTATUS[@]}")
    if [ "${pipe_status[0]}" != "0" ] || [ ! -s ${BAM}.tmp ]; then
        echo "    ⚠ mapping failed (status: ${pipe_status[*]})"
        rm -f ${BAM}.tmp
        continue
    fi
    
    mv ${BAM}.tmp $BAM
    samtools index $BAM
    
    # 简短报告
    n_mapped=$(samtools view -c -F 4 $BAM)
    echo "    mapped reads: $n_mapped"
done < $SAMPLES_TSV

echo ""
echo "BAMs ready: $(ls -1 ${BAM_DIR}/*.bam 2>/dev/null | wc -l) / 26"

# -------------------- Step 3: inStrain profile --------------------
echo ""
echo "[$(date +%H:%M:%S)] Running inStrain profile..."

cd $BAM_DIR
for filename in *.bam; do
    base=$(basename $filename .bam)
    out=../${PROFILE_DIR}/${base}_profile
    
    if [ -d "$out" ] && [ -f "${out}/output/${base}_profile_genome_info.tsv" ]; then
        echo "  [done] $base"
        continue
    fi
    
    echo "  [run]  $base"
    inStrain profile \
        ${base}.bam $REF \
        -s $STB \
        --gene_file $GENES \
        -p $THREADS \
        -o $out \
        > ../${LOG_DIR}/instrain_${base}.log 2>&1 || \
        echo "    ⚠ failed (see ${LOG_DIR}/instrain_${base}.log)"
done
cd ..

# -------------------- Step 4: 汇总（per-sample × per-genome） --------------------
echo ""
echo "[$(date +%H:%M:%S)] 汇总每个 sample × 每个 genome 的统计..."

OUT=results/instrain_per_sample_per_genome.tsv
header_done=0
> $OUT

while IFS=$'\t' read -r sample clade rest; do
    [[ $sample == "sample_id" ]] && continue
    
    info="${PROFILE_DIR}/${sample}_profile/output/${sample}_profile_genome_info.tsv"
    [ -f "$info" ] || continue
    
    if [ $header_done -eq 0 ]; then
        head -1 $info | awk -v OFS='\t' '{print "sample","clade",$0}' > $OUT
        header_done=1
    fi
    
    awk -v s=$sample -v c=$clade -v OFS='\t' 'NR>1 {print s,c,$0}' $info >> $OUT
done < $SAMPLES_TSV

echo "  → $OUT  ($(wc -l < $OUT) 行)"

# -------------------- 总结 --------------------
echo ""
echo "=========================================="
echo " 完成"
echo "=========================================="
echo ""
echo "输出："
echo "  ${BAM_DIR}/                              26 个 sample 的 BAM (map to all_AOA_contig)"
echo "  ${PROFILE_DIR}/<sample>_profile/         每个 sample 的 inStrain profile"
echo "  results/instrain_per_sample_per_genome.tsv   汇总 (sample × genome)"
echo ""
echo "下一步可选："
echo "  - inStrain compare -i ${PROFILE_DIR}/*_profile -o instrain_compare -p $THREADS"
echo "    (做 pairwise popANI / conANI)"
