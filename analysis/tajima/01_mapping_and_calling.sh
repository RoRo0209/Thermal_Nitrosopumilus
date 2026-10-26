#!/bin/bash
: "${THERMAL_PROJECT_ROOT:?Set THERMAL_PROJECT_ROOT to the analysis workspace}"
cd "${THERMAL_PROJECT_ROOT}/tajima"

# ==============================================================================
# 01_mapping_and_calling.sh (v3 - 更稳健)
# 改动：
#   1. 不用 subshell，用显式的 PIPESTATUS 检查，确保错误处理到位
#   2. 默认线程数降到 8，避免把共享服务器打爆
#   3. 每个样本跑之前检查服务器负载
#   4. Breadth < 20% 自动跳过 variant calling
#   5. freebayes 加 timeout
# ==============================================================================

set -uo pipefail

# -------------------- 配置 --------------------
MANIFEST="sample_manifest.tsv"
HIGHT_REF="refs/HighT_representative.fna"
LOWT_REF="refs/LowT_representative.fna"
THREADS=8                        # 降到 8（之前是 16）
MIN_MAPQ=20
MIN_ALT_FRAC=0.05
MIN_ALT_COUNT=2
MIN_COVERAGE=3
MIN_BREADTH_FOR_VC=0.2           # breadth < 20% 不做 variant calling
MAX_LOAD_MULT=4                  # load average > CPU核数 * 4 就等待
FREEBAYES_TIMEOUT=3600            # freebayes 最多跑 1 小时

mkdir -p bam vcf logs coverage

# -------------------- 失败记录 --------------------
FAILED_LOG="failed_samples.tsv"
if [ ! -f $FAILED_LOG ]; then
    echo -e "sample_id\tclade\tstep\treason\ttimestamp" > $FAILED_LOG
fi

record_failure() {
    local sample=$1 clade=$2 step=$3 reason=$4
    echo -e "${sample}\t${clade}\t${step}\t${reason}\t$(date +%F_%T)" >> $FAILED_LOG
    echo "  ✗ FAILED: $sample ($clade) at $step - $reason"
}

# -------------------- 负载检查 --------------------
wait_for_load() {
    local ncpu=$(nproc)
    local max_load=$((ncpu * MAX_LOAD_MULT))
    while true; do
        local load=$(awk '{print $1}' /proc/loadavg | awk -F. '{print $1}')
        if [ $load -lt $max_load ]; then
            return 0
        fi
        echo "  ⏸ Load average = $load (max allowed: $max_load), waiting 60s..."
        sleep 60
    done
}

# -------------------- 预检查 reads --------------------
check_paired_reads() {
    local r1=$1 r2=$2
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        echo "MISSING_FILE"; return 1
    fi
    local s1=$(stat -c%s "$r1")
    local s2=$(stat -c%s "$r2")
    local ratio=$(awk -v a=$s1 -v b=$s2 'BEGIN{print (a>b)? a/b : b/a}')
    local big=$(awk -v r=$ratio 'BEGIN{print (r>1.5)? 1 : 0}')
    if [ "$big" = "1" ]; then
        echo "SIZE_MISMATCH"; return 1
    fi
    echo "OK"; return 0
}

# -------------------- 建索引 --------------------
for REF in $HIGHT_REF $LOWT_REF; do
    if [ ! -f ${REF%.fna}.1.bt2 ] && [ ! -f ${REF%.fna}.1.bt2l ]; then
        echo "[$(date +%H:%M:%S)] Building bowtie2 index for $REF"
        bowtie2-build $REF ${REF%.fna}
    fi
    if [ ! -f ${REF}.fai ]; then
        samtools faidx $REF
    fi
done

# -------------------- 主循环 --------------------
tail -n +2 $MANIFEST | while IFS=$'\t' read -r sample r1 r2 highT_frac lowT_frac clade; do
    
    [[ $sample == \#* ]] && continue
    [[ $clade == "mixed" || $clade == "exclude" || $clade == "no_signal" ]] && continue
    
    if [[ $clade == "HighT" ]]; then
        REF=$HIGHT_REF; TAG="HighT"
    elif [[ $clade == "LowT" ]]; then
        REF=$LOWT_REF; TAG="LowT"
    else
        continue
    fi
    
    REF_BASE=${REF%.fna}
    BAM="bam/${sample}_${TAG}.bam"
    VCF="vcf/${sample}_${TAG}.vcf.gz"
    
    # 断点续跑
    if [ -f $VCF ] && [ -f coverage/${sample}_${TAG}.cov ]; then
        continue
    fi
    
    echo ""
    echo "[$(date +%H:%M:%S)] ===== $sample ($TAG) ====="
    
    # 等待服务器负载降低
    wait_for_load
    
    # -------- Step 0: 预检查 reads --------
    check_result=$(check_paired_reads "$r1" "$r2")
    if [ "$check_result" != "OK" ]; then
        record_failure $sample $TAG "pre_check" "$check_result"
        continue
    fi
    
    # -------- Step 1: Mapping --------
    if [ ! -f $BAM ]; then
        echo "  [1/3] Mapping..."
        
        # 用显式的错误检查，不用 subshell
        bowtie2 -x $REF_BASE \
            -1 $r1 -2 $r2 \
            --very-sensitive --no-unal \
            -p $THREADS \
            2> logs/${sample}_${TAG}_bowtie.log | \
        samtools view -b -q $MIN_MAPQ -F 4 -@ 2 - | \
        samtools sort -@ 2 -o $BAM.tmp -
        
        # 检查管道里每一步的退出状态
        pipe_status=("${PIPESTATUS[@]}")
        failed=0
        for s in "${pipe_status[@]}"; do
            if [ "$s" != "0" ]; then
                failed=1
                break
            fi
        done
        
        if [ $failed -eq 1 ] || [ ! -s $BAM.tmp ]; then
            record_failure $sample $TAG "mapping" "bowtie2/samtools exit=${pipe_status[*]}"
            rm -f $BAM.tmp
            continue
        fi
        
        mv $BAM.tmp $BAM
        samtools index $BAM
    fi
    
    # -------- Step 2: Coverage check --------
    if [ ! -f coverage/${sample}_${TAG}.cov ]; then
        echo "  [2/3] Coverage check..."
        depth=$(samtools depth -a $BAM 2>/dev/null | awk '{sum+=$3; n++} END{if(n>0) print sum/n; else print 0}')
        breadth=$(samtools depth -a $BAM 2>/dev/null | awk -v m=$MIN_COVERAGE '$3>=m{n++} END{if(NR>0) print n/NR; else print 0}')
        echo -e "${sample}\t${TAG}\t${depth}\t${breadth}" > coverage/${sample}_${TAG}.cov
        echo "    depth=${depth}, breadth@${MIN_COVERAGE}x=${breadth}"
    else
        depth=$(awk '{print $3}' coverage/${sample}_${TAG}.cov)
        breadth=$(awk '{print $4}' coverage/${sample}_${TAG}.cov)
    fi
    
    # Breadth 太低直接跳
    breadth_ok=$(awk -v b=$breadth -v t=$MIN_BREADTH_FOR_VC 'BEGIN{print (b >= t)? 1 : 0}')
    if [ "$breadth_ok" = "0" ]; then
        record_failure $sample $TAG "low_breadth" "breadth=${breadth}<${MIN_BREADTH_FOR_VC}"
        continue
    fi
    
    # -------- Step 3: Variant calling --------
    if [ ! -f $VCF ]; then
        echo "  [3/3] Variant calling..."
        
        timeout $FREEBAYES_TIMEOUT freebayes -f $REF \
            --ploidy 1 \
            --pooled-continuous \
            --min-alternate-fraction $MIN_ALT_FRAC \
            --min-alternate-count $MIN_ALT_COUNT \
            --min-coverage $MIN_COVERAGE \
            --report-monomorphic \
            $BAM 2> logs/${sample}_${TAG}_freebayes.log | \
        bcftools view -Oz -o ${VCF}.tmp 2>> logs/${sample}_${TAG}_freebayes.log
        
        pipe_status=("${PIPESTATUS[@]}")
        fb_status=${pipe_status[0]}
        bcf_status=${pipe_status[1]}
        
        # 检查 freebayes 和 bcftools 是否都成功
        if [ "$fb_status" != "0" ]; then
            if [ "$fb_status" = "124" ]; then
                record_failure $sample $TAG "variant_calling" "freebayes_timeout"
            else
                record_failure $sample $TAG "variant_calling" "freebayes_exit=${fb_status}"
            fi
            rm -f ${VCF}.tmp
            continue
        fi
        
        if [ "$bcf_status" != "0" ]; then
            record_failure $sample $TAG "variant_calling" "bcftools_exit=${bcf_status}"
            rm -f ${VCF}.tmp
            continue
        fi
        
        if [ ! -s ${VCF}.tmp ]; then
            record_failure $sample $TAG "variant_calling" "empty_vcf"
            rm -f ${VCF}.tmp
            continue
        fi
        
        mv ${VCF}.tmp $VCF
        tabix -p vcf $VCF
    fi
    
    echo "  ✓ DONE"
done

# -------------------- 汇总 --------------------
echo ""
echo "=========================================="
echo " Pipeline Summary"
echo "=========================================="

echo -e "sample\tclade\tdepth\tbreadth" > coverage_summary.tsv
cat coverage/*.cov 2>/dev/null >> coverage_summary.tsv

n_done=$(ls vcf/*.vcf.gz 2>/dev/null | wc -l)
n_failed=$(tail -n +2 $FAILED_LOG 2>/dev/null | wc -l)
n_total=$(awk -F'\t' 'NR>1 && ($6=="HighT" || $6=="LowT")' $MANIFEST | wc -l)

echo "Total mapping tasks: $n_total"
echo "Completed:           $n_done"
echo "Failed:              $n_failed"

if [ $n_failed -gt 0 ]; then
    echo ""
    echo "Failed samples (see failed_samples.tsv):"
    tail -n +2 $FAILED_LOG | column -t
fi

echo ""
echo "Top 20 samples by breadth:"
sort -k4 -rn coverage_summary.tsv 2>/dev/null | head -20 | column -t
