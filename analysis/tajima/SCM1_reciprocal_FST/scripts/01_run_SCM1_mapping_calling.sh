#!/usr/bin/env bash
: "${THERMAL_PROJECT_ROOT:?Set THERMAL_PROJECT_ROOT to the analysis workspace}"
cd "${THERMAL_PROJECT_ROOT}/tajima/SCM1_reciprocal_FST"

# Reciprocal validation: map the same 26 FST samples to the HighT representative
# genome SCM1 and generate a shared-coordinate joint VCF.

set -euo pipefail

BASE_DIR="${THERMAL_PROJECT_ROOT}/tajima/SCM1_reciprocal_FST"
TAJIMA_DIR="${THERMAL_PROJECT_ROOT}/tajima"
cd "${BASE_DIR}"

MANIFEST="${TAJIMA_DIR}/sample_manifest.tsv"
SAMPLE_GROUPS="${BASE_DIR}/sample_groups_26.tsv"
REF="${TAJIMA_DIR}/refs/HighT_representative.fna"
REF_PREFIX="${TAJIMA_DIR}/refs/HighT_representative"

THREADS="${THREADS:-8}"
MIN_MAPQ="${MIN_MAPQ:-20}"
MIN_ALT_FRAC="${MIN_ALT_FRAC:-0.05}"
MIN_ALT_COUNT="${MIN_ALT_COUNT:-2}"
MIN_COVERAGE="${MIN_COVERAGE:-3}"
QC_MIN_DEPTH="${QC_MIN_DEPTH:-5}"
QC_MIN_BREADTH="${QC_MIN_BREADTH:-0.3}"

mkdir -p bam vcf logs coverage results tmp

echo "[INFO] SCM1 reciprocal mapping started: $(date)"
echo "[INFO] BASE_DIR=${BASE_DIR}"
echo "[INFO] REF=${REF}"

for tool in bowtie2 samtools freebayes bcftools tabix awk; do
    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: ${tool}" >&2
        exit 2
    fi
done

if [[ ! -s "${REF}" ]]; then
    echo "[ERROR] SCM1 reference fasta not found: ${REF}" >&2
    exit 2
fi

if [[ ! -s "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi

if [[ ! -s "${REF_PREFIX}.1.bt2" && ! -s "${REF_PREFIX}.1.bt2l" ]]; then
    bowtie2-build "${REF}" "${REF_PREFIX}"
fi

printf "sample\tgroup\tr1\tr2\tbam\tvcf\tmean_depth\tbreadth_ge%d\tpasses_original_thresholds\n" "${MIN_COVERAGE}" > results/SCM1_mapping_QC.tsv

tail -n +2 "${SAMPLE_GROUPS}" | while IFS=$'\t' read -r sample group; do
    r1="$(awk -F'\t' -v s="${sample}" 'NR>1 && $1==s {print $2; exit}' "${MANIFEST}")"
    r2="$(awk -F'\t' -v s="${sample}" 'NR>1 && $1==s {print $3; exit}' "${MANIFEST}")"
    if [[ -z "${r1}" || -z "${r2}" ]]; then
        echo "[ERROR] Reads not found in manifest for ${sample}" >&2
        exit 3
    fi
    if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
        echo "[ERROR] Read files not accessible for ${sample}: ${r1} ${r2}" >&2
        exit 3
    fi

    bam="bam/${sample}.SCM1.sorted.bam"
    vcf="vcf/${sample}.SCM1.vcf.gz"

    echo "[INFO] Mapping ${sample} (${group})"
    if [[ ! -s "${bam}" ]]; then
        bowtie2 -x "${REF_PREFIX}" \
            -1 "${r1}" -2 "${r2}" \
            --very-sensitive --no-unal \
            -p "${THREADS}" 2> "logs/${sample}.bowtie2.log" |
        samtools view -b -q "${MIN_MAPQ}" -F 4 -@ 2 - |
        samtools sort -@ 2 -o "${bam}" -
        samtools index "${bam}"
    fi

    depth="$(samtools depth -a "${bam}" | awk '{sum+=$3; n++} END{if(n>0) print sum/n; else print 0}')"
    breadth="$(samtools depth -a "${bam}" | awk -v m="${MIN_COVERAGE}" '$3>=m{n++} END{if(NR>0) print n/NR; else print 0}')"
    pass="$(awk -v d="${depth}" -v b="${breadth}" -v md="${QC_MIN_DEPTH}" -v mb="${QC_MIN_BREADTH}" 'BEGIN{print (d>=md && b>=mb) ? "yes" : "no"}')"

    echo "[INFO] ${sample}: depth=${depth}, breadth=${breadth}, pass=${pass}"

    if [[ ! -s "${vcf}" ]]; then
        freebayes -f "${REF}" \
            --ploidy 1 --pooled-continuous \
            --min-alternate-fraction "${MIN_ALT_FRAC}" \
            --min-alternate-count "${MIN_ALT_COUNT}" \
            --min-coverage "${MIN_COVERAGE}" \
            --report-monomorphic \
            "${bam}" 2> "logs/${sample}.freebayes.log" |
        bcftools view -Oz -o "${vcf}.tmp"

        echo "${sample}" > "tmp/${sample}.name.txt"
        bcftools reheader --samples "tmp/${sample}.name.txt" "${vcf}.tmp" -o "${vcf}"
        tabix -f -p vcf "${vcf}"
        rm -f "${vcf}.tmp" "tmp/${sample}.name.txt"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${sample}" "${group}" "${r1}" "${r2}" "${bam}" "${vcf}" "${depth}" "${breadth}" "${pass}" \
        >> results/SCM1_mapping_QC.tsv
done

echo "[INFO] Merging per-sample VCF files"
vcf_list="tmp/scm1_vcfs.list"
tail -n +2 "${SAMPLE_GROUPS}" | awk -v base="${BASE_DIR}" '{print base"/vcf/"$1".SCM1.vcf.gz"}' > "${vcf_list}"

bcftools merge --file-list "${vcf_list}" \
    -m none --missing-to-ref \
    -Oz -o results/SCM1_joint_26samples.vcf.gz
tabix -f -p vcf results/SCM1_joint_26samples.vcf.gz

echo "[INFO] Computing amoABC per-sample depth"
printf "sample\tgroup\tgene\tcontig\tstart\tend\tmean_depth\tbreadth_ge%d\n" "${MIN_COVERAGE}" > results/SCM1_amoABC_depth_by_sample.tsv
tail -n +2 SCM1_amoABC_loci.tsv | while IFS=$'\t' read -r gene locus gff_seqid contig start end strand source note; do
    tail -n +2 "${SAMPLE_GROUPS}" | while IFS=$'\t' read -r sample group; do
        bam="bam/${sample}.SCM1.sorted.bam"
        region="${contig}:${start}-${end}"
        mean_depth="$(samtools depth -a -r "${region}" "${bam}" | awk '{sum+=$3; n++} END{if(n>0) print sum/n; else print 0}')"
        breadth="$(samtools depth -a -r "${region}" "${bam}" | awk -v m="${MIN_COVERAGE}" '$3>=m{n++} END{if(NR>0) print n/NR; else print 0}')"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${sample}" "${group}" "${gene}" "${contig}" "${start}" "${end}" "${mean_depth}" "${breadth}" \
            >> results/SCM1_amoABC_depth_by_sample.tsv
    done
done

echo "[INFO] Finished SCM1 reciprocal mapping/calling: $(date)"
echo "[INFO] Next: python scripts/02_analyze_SCM1_reciprocal.py"
