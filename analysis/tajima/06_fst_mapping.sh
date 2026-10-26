#!/bin/bash
: "${THERMAL_PROJECT_ROOT:?Set THERMAL_PROJECT_ROOT to the analysis workspace}"
cd "${THERMAL_PROJECT_ROOT}/tajima"

# ==============================================================================
# 06_fst_analysis.sh
# 功能：把 HighT 的达标样本重 map 到 LowT reference，用于算 Fst
# ==============================================================================

set -uo pipefail

# -------------------- 配置 --------------------
# 只处理 HighT 样本（LowT 已经 map 到 LowT ref 了）
MANIFEST="sample_manifest.tsv"
SHARED_REF="refs/LowT_representative.fna"   # 用 LowT 作共同 reference
THREADS=8
MIN_MAPQ=20
MIN_ALT_FRAC=0.05
MIN_ALT_COUNT=2
MIN_COVERAGE=3

mkdir -p bam_shared vcf_shared logs_shared coverage_shared

# 建索引（如果还没有）
if [ ! -f ${SHARED_REF%.fna}.1.bt2 ] && [ ! -f ${SHARED_REF%.fna}.1.bt2l ]; then
    bowtie2-build $SHARED_REF ${SHARED_REF%.fna}
fi

# 用 HighT 达标样本清单
if [ ! -f HighT_samples_passing.txt ]; then
    echo "错误：需要先跑 02_merge_populations.sh 生成 HighT_samples_passing.txt"
    exit 1
fi

# -------------------- 对每个 HighT 样本重 map --------------------
echo "[$(date +%H:%M:%S)] 重新 mapping HighT 样本到 LowT reference..."

while read sample; do
    # 从 manifest 找 reads 路径
    read r1 r2 <<< $(awk -F'\t' -v s=$sample 'NR>1 && $1==s {print $2" "$3}' $MANIFEST)
    
    if [ -z "$r1" ] || [ -z "$r2" ]; then
        echo "  ⚠ $sample: reads 路径没找到，跳过"
        continue
    fi
    
    BAM="bam_shared/${sample}_sharedref.bam"
    VCF="vcf_shared/${sample}_sharedref.vcf.gz"
    
    if [ -f $VCF ]; then
        echo "  [SKIP] $sample 已完成"
        continue
    fi
    
    echo "  [mapping] $sample"
    bowtie2 -x ${SHARED_REF%.fna} \
        -1 $r1 -2 $r2 \
        --very-sensitive --no-unal \
        -p $THREADS 2> logs_shared/${sample}_bowtie.log | \
    samtools view -b -q $MIN_MAPQ -F 4 -@ 2 - | \
    samtools sort -@ 2 -o $BAM -
    samtools index $BAM
    
    # Coverage check
    depth=$(samtools depth -a $BAM 2>/dev/null | awk '{sum+=$3; n++} END{if(n>0) print sum/n; else print 0}')
    breadth=$(samtools depth -a $BAM 2>/dev/null | awk -v m=$MIN_COVERAGE '$3>=m{n++} END{if(NR>0) print n/NR; else print 0}')
    echo -e "${sample}\tHighT_on_LowT\t${depth}\t${breadth}" > coverage_shared/${sample}.cov
    echo "    depth=${depth}, breadth=${breadth}"
    
    # Variant calling
    echo "  [variant calling] $sample"
    timeout 3600 freebayes -f $SHARED_REF \
        --ploidy 1 --pooled-continuous \
        --min-alternate-fraction $MIN_ALT_FRAC \
        --min-alternate-count $MIN_ALT_COUNT \
        --min-coverage $MIN_COVERAGE \
        --report-monomorphic \
        $BAM 2> logs_shared/${sample}_freebayes.log | \
    bcftools view -Oz -o ${VCF}.tmp
    
    if [ -s ${VCF}.tmp ]; then
        # 改 sample name
        echo "$sample" > /tmp/newname_$$.txt
        bcftools reheader --samples /tmp/newname_$$.txt ${VCF}.tmp -o $VCF
        tabix -f -p vcf $VCF
        rm -f ${VCF}.tmp /tmp/newname_$$.txt
    else
        rm -f ${VCF}.tmp
        echo "  ⚠ $sample: VCF empty, skipping"
    fi
done < HighT_samples_passing.txt

# -------------------- 合并 HighT VCF --------------------
echo ""
echo "[$(date +%H:%M:%S)] 合并 HighT on shared reference..."

HIGHT_SHARED_VCFS=""
while read sample; do
    vcf="vcf_shared/${sample}_sharedref.vcf.gz"
    if [ -f $vcf ]; then
        HIGHT_SHARED_VCFS="$HIGHT_SHARED_VCFS $vcf"
    fi
done < HighT_samples_passing.txt

bcftools merge $HIGHT_SHARED_VCFS \
    -m none --missing-to-ref \
    -Oz -o vcf_shared/HighT_on_LowTref.vcf.gz
tabix -p vcf vcf_shared/HighT_on_LowTref.vcf.gz

echo "  ✓ vcf_shared/HighT_on_LowTref.vcf.gz"

# -------------------- 合并 HighT + LowT 为 joint VCF --------------------
echo ""
echo "[$(date +%H:%M:%S)] 创建 joint VCF (HighT + LowT on same reference)..."

bcftools merge vcf_shared/HighT_on_LowTref.vcf.gz vcf/LowT_population.vcf.gz \
    -m none --missing-to-ref \
    -Oz -o vcf_shared/HighT_LowT_joint.vcf.gz
tabix -p vcf vcf_shared/HighT_LowT_joint.vcf.gz

echo "  ✓ vcf_shared/HighT_LowT_joint.vcf.gz"
echo ""
echo "样本组成:"
bcftools query -l vcf_shared/HighT_LowT_joint.vcf.gz

echo ""
echo "=========================================="
echo " 完成！下一步运行 python 07_fst_analysis.py"
echo "=========================================="
