#!/usr/bin/env bash
: "${THERMAL_PROJECT_ROOT:?Set THERMAL_PROJECT_ROOT to the analysis workspace}"
cd "${THERMAL_PROJECT_ROOT}/tajima"

set -euo pipefail

# Call variants only in the SPOT01 C2 / amoABC region from workstation BAMs.
# Output is intentionally small and can be copied back to the local project.
#
# Usage:
#   bash workstation_amoabc_call_and_export.sh \
#     /path/to/bam_shared:/path/to/bam \
#     refs/LowT_representative.fna \
#     sample_manifest.tsv \
#     amoabc_workstation_out
#
# Assumptions:
#   - BAMs are mapped to SPOT01 / LowT reference.
#   - BAM filenames contain sample IDs from sample_manifest.tsv.
#   - BAM_DIR can be one directory, or multiple directories separated by ":" or ",".
#   - Optional: set HIGHT_BAM_DIRS and LOWT_BAM_DIRS to search different
#     directories by clade. This is useful when HighT shared-reference BAMs
#     are in bam_shared/ and LowT SPOT01 BAMs are in bam/.
#   - bcftools, samtools and tabix are available on the workstation.

BAM_DIRS=${1:?BAM directory or directory list required}
REF=${2:?Reference FASTA required}
MANIFEST=${3:?sample_manifest.tsv required}
OUTDIR=${4:-amoabc_workstation_out}

REGION="${THERMAL_AMO_REGION:?Set THERMAL_AMO_REGION to contig:start-end on the matched reference}"
MIN_MAPQ=${MIN_MAPQ:-20}
MIN_BASEQ=${MIN_BASEQ:-20}
THREADS=${THREADS:-8}

mkdir -p "$OUTDIR"

if [ ! -f "${REF}.fai" ]; then
  samtools faidx "$REF"
fi

echo "[1/5] Reading qualified HighT/LowT samples from $MANIFEST"
awk -F'\t' 'NR>1 && ($6=="HighT" || $6=="LowT") {print $1 "\t" $6}' "$MANIFEST" \
  > "$OUTDIR/sample_groups.tsv"

echo "[2/5] Matching BAM files in $BAM_DIRS"
: > "$OUTDIR/bams.txt"
: > "$OUTDIR/missing_bams.tsv"

IFS=':,' read -r -a BAM_DIR_ARRAY <<< "$BAM_DIRS"
IFS=':,' read -r -a HIGHT_BAM_DIR_ARRAY <<< "${HIGHT_BAM_DIRS:-$BAM_DIRS}"
IFS=':,' read -r -a LOWT_BAM_DIR_ARRAY <<< "${LOWT_BAM_DIRS:-$BAM_DIRS}"

while IFS=$'\t' read -r sample clade; do
  bam=""
  if [ "$clade" = "HighT" ]; then
    SEARCH_DIRS=("${HIGHT_BAM_DIR_ARRAY[@]}")
  elif [ "$clade" = "LowT" ]; then
    SEARCH_DIRS=("${LOWT_BAM_DIR_ARRAY[@]}")
  else
    SEARCH_DIRS=("${BAM_DIR_ARRAY[@]}")
  fi
  for bam_dir in "${SEARCH_DIRS[@]}"; do
    [ -d "$bam_dir" ] || continue
    hit=$(find "$bam_dir" -maxdepth 1 -type f -name "${sample}*.bam" | head -n 1)
    if [ -n "$hit" ]; then
      bam="$hit"
      break
    fi
  done
  if [ -z "$bam" ]; then
    echo -e "${sample}\t${clade}" >> "$OUTDIR/missing_bams.tsv"
    continue
  fi
  if [ ! -f "${bam}.bai" ] && [ ! -f "${bam%.bam}.bai" ]; then
    samtools index "$bam"
  fi
  echo "$bam" >> "$OUTDIR/bams.txt"
done < "$OUTDIR/sample_groups.tsv"

N_BAM=$(wc -l < "$OUTDIR/bams.txt" | tr -d ' ')
echo "  Found $N_BAM BAMs"

if [ "$N_BAM" -lt 4 ]; then
  echo "Too few BAMs found. Check BAM_DIR and sample names." >&2
  exit 1
fi

echo "[3/5] Exporting per-sample depth over $REGION"
samtools depth -a -q "$MIN_MAPQ" -Q "$MIN_BASEQ" -r "$REGION" -f "$OUTDIR/bams.txt" \
  > "$OUTDIR/amoabc.depth.tsv"

echo "[4/5] Calling small joint VCF over $REGION"
bcftools mpileup \
  --threads "$THREADS" \
  -f "$REF" \
  -r "$REGION" \
  -q "$MIN_MAPQ" \
  -Q "$MIN_BASEQ" \
  -a FORMAT/AD,FORMAT/DP \
  -Ou \
  -b "$OUTDIR/bams.txt" \
| bcftools call \
  --threads "$THREADS" \
  -m -A -v \
  -Oz \
  -o "$OUTDIR/amoabc.raw.vcf.gz"

tabix -f -p vcf "$OUTDIR/amoabc.raw.vcf.gz"

echo "[5/5] Normalizing multiallelic sites"
bcftools norm \
  --threads "$THREADS" \
  -m-any \
  -f "$REF" \
  -Oz \
  -o "$OUTDIR/amoabc.norm.vcf.gz" \
  "$OUTDIR/amoabc.raw.vcf.gz"

tabix -f -p vcf "$OUTDIR/amoabc.norm.vcf.gz"

bcftools query -l "$OUTDIR/amoabc.norm.vcf.gz" > "$OUTDIR/vcf_samples.txt"
bcftools view -H "$OUTDIR/amoabc.norm.vcf.gz" | wc -l > "$OUTDIR/n_variant_sites.txt"

echo "Done."
echo "Copy these small files back:"
echo "  $OUTDIR/amoabc.norm.vcf.gz"
echo "  $OUTDIR/amoabc.norm.vcf.gz.tbi"
echo "  $OUTDIR/sample_groups.tsv"
echo "  $OUTDIR/vcf_samples.txt"
echo "  $OUTDIR/n_variant_sites.txt"
