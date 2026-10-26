#!/bin/bash
: "${THERMAL_PROJECT_ROOT:?Set THERMAL_PROJECT_ROOT to the analysis workspace}"
cd "${THERMAL_PROJECT_ROOT}/tajima"

# ==============================================================================
# 02_merge_populations.sh (v2)
# 改进：
#   1. 先给每个 per-sample VCF 加正确的 sample name (从文件名提取)
#   2. 再 merge，不会冲突
#   3. 检查最终合并 VCF 的内容
# ==============================================================================

set -uo pipefail

# -------------------- 配置 --------------------
MIN_DEPTH=5
MIN_BREADTH=0.4

mkdir -p vcf_renamed

# -------------------- Step 1: 按 clade 筛选达标样本 --------------------
echo "[$(date +%H:%M:%S)] 筛选达标样本..."

awk -v d=$MIN_DEPTH -v b=$MIN_BREADTH \
    'NR>1 && $3>=d && $4>=b && $2=="HighT" {print $1}' \
    coverage_summary.tsv > HighT_samples_passing.txt

awk -v d=$MIN_DEPTH -v b=$MIN_BREADTH \
    'NR>1 && $3>=d && $4>=b && $2=="LowT" {print $1}' \
    coverage_summary.tsv > LowT_samples_passing.txt

N_HIGHT=$(wc -l < HighT_samples_passing.txt)
N_LOWT=$(wc -l < LowT_samples_passing.txt)

echo "  HighT population: $N_HIGHT samples pass"
echo "  LowT population:  $N_LOWT samples pass"

if [ $N_HIGHT -lt 3 ] || [ $N_LOWT -lt 3 ]; then
    echo "⚠ 警告：某个 clade 的合格样本数 < 3"
    echo "     你可以降低 MIN_DEPTH 或 MIN_BREADTH 再跑"
    echo "     当前: MIN_DEPTH=$MIN_DEPTH, MIN_BREADTH=$MIN_BREADTH"
    if [ $N_HIGHT -lt 1 ] || [ $N_LOWT -lt 1 ]; then
        echo "致命：完全没样本，退出"
        exit 1
    fi
fi

# -------------------- Step 2: 给每个 VCF 重命名 sample name --------------------
echo ""
echo "[$(date +%H:%M:%S)] 给 VCF 改 sample name..."

rename_vcf_sample() {
    local in_vcf=$1
    local new_name=$2
    local out_vcf=$3
    
    # 生成 sample name 文件
    echo "$new_name" > /tmp/newname_$$.txt
    
    # 用 bcftools reheader 改名
    bcftools reheader --samples /tmp/newname_$$.txt "$in_vcf" -o "$out_vcf"
    tabix -f -p vcf "$out_vcf"
    
    rm -f /tmp/newname_$$.txt
}

# 处理 HighT
HIGHT_RENAMED_VCFS=""
while read sample; do
    in_vcf="vcf/${sample}_HighT.vcf.gz"
    out_vcf="vcf_renamed/${sample}_HighT.vcf.gz"
    
    if [ ! -f "$in_vcf" ]; then
        echo "  ⚠ 缺少 $in_vcf，跳过"
        continue
    fi
    
    if [ ! -f "$out_vcf" ]; then
        rename_vcf_sample "$in_vcf" "$sample" "$out_vcf"
        echo "  ✓ $sample (HighT)"
    fi
    HIGHT_RENAMED_VCFS="$HIGHT_RENAMED_VCFS $out_vcf"
done < HighT_samples_passing.txt

# 处理 LowT
LOWT_RENAMED_VCFS=""
while read sample; do
    in_vcf="vcf/${sample}_LowT.vcf.gz"
    out_vcf="vcf_renamed/${sample}_LowT.vcf.gz"
    
    if [ ! -f "$in_vcf" ]; then
        echo "  ⚠ 缺少 $in_vcf，跳过"
        continue
    fi
    
    if [ ! -f "$out_vcf" ]; then
        rename_vcf_sample "$in_vcf" "$sample" "$out_vcf"
        echo "  ✓ $sample (LowT)"
    fi
    LOWT_RENAMED_VCFS="$LOWT_RENAMED_VCFS $out_vcf"
done < LowT_samples_passing.txt

# -------------------- Step 3: 合并 --------------------
echo ""
echo "[$(date +%H:%M:%S)] 合并 HighT population VCF..."

if [ -z "$HIGHT_RENAMED_VCFS" ]; then
    echo "  ⚠ 没有 HighT VCF 可合并"
else
    bcftools merge $HIGHT_RENAMED_VCFS \
        -m none \
        --missing-to-ref \
        -Oz -o vcf/HighT_population.vcf.gz
    tabix -f -p vcf vcf/HighT_population.vcf.gz
    echo "  ✓ vcf/HighT_population.vcf.gz"
fi

echo ""
echo "[$(date +%H:%M:%S)] 合并 LowT population VCF..."

if [ -z "$LOWT_RENAMED_VCFS" ]; then
    echo "  ⚠ 没有 LowT VCF 可合并"
else
    bcftools merge $LOWT_RENAMED_VCFS \
        -m none \
        --missing-to-ref \
        -Oz -o vcf/LowT_population.vcf.gz
    tabix -f -p vcf vcf/LowT_population.vcf.gz
    echo "  ✓ vcf/LowT_population.vcf.gz"
fi

# -------------------- Step 4: 生成 pixy populations 文件 --------------------
awk 'BEGIN{OFS="\t"} {print $1, "HighT"}' HighT_samples_passing.txt > pops_HighT.txt
awk 'BEGIN{OFS="\t"} {print $1, "LowT"}' LowT_samples_passing.txt > pops_LowT.txt

# -------------------- Step 5: 汇总 --------------------
echo ""
echo "=========================================="
echo " 完成！"
echo "=========================================="

if [ -f vcf/HighT_population.vcf.gz ]; then
    n_snv_h=$(bcftools view -H vcf/HighT_population.vcf.gz 2>/dev/null | wc -l)
    echo ""
    echo "HighT population VCF: vcf/HighT_population.vcf.gz"
    echo "  样本数: $N_HIGHT"
    echo "  位点数: $n_snv_h"
    # 显示 sample names 确认没重复
    echo "  sample names:"
    bcftools query -l vcf/HighT_population.vcf.gz | sed 's/^/    /'
fi

if [ -f vcf/LowT_population.vcf.gz ]; then
    n_snv_l=$(bcftools view -H vcf/LowT_population.vcf.gz 2>/dev/null | wc -l)
    echo ""
    echo "LowT population VCF:  vcf/LowT_population.vcf.gz"
    echo "  样本数: $N_LOWT"
    echo "  位点数: $n_snv_l"
    echo "  sample names:"
    bcftools query -l vcf/LowT_population.vcf.gz | sed 's/^/    /'
fi
