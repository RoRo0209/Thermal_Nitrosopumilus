#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 INPUT_05_DIR OUTPUT_07_DIR [THREADS]" >&2
  exit 2
fi

INPUT_DIR="$(realpath "$1")"
OUTPUT_DIR="$2"
THREADS="${3:-4}"
RESUME="${RESUME:-1}"

PROTEIN_DIR="$INPUT_DIR/qc_proteins"
CDS_DIR="$INPUT_DIR/qc_cds"
mkdir -p "$OUTPUT_DIR/protein_aligned" "$OUTPUT_DIR/protein_trimmed" \
         "$OUTPUT_DIR/codon_full" "$OUTPUT_DIR/codon_aligned" "$OUTPUT_DIR/logs"

for tool in mafft pal2nal.pl trimal; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Required executable not found: $tool" >&2
    exit 1
  fi
done

printf 'file_stem\tstatus\n' > "$OUTPUT_DIR/alignment_status.tsv"
processed=0
skipped=0

shopt -s nullglob
for protein in "$PROTEIN_DIR"/*.faa; do
  stem="$(basename "$protein" .faa)"
  cds="$CDS_DIR/$stem.fna"
  aligned="$OUTPUT_DIR/protein_aligned/$stem.aln.faa"
  trimmed="$OUTPUT_DIR/protein_trimmed/$stem.trim.faa"
  codon_full="$OUTPUT_DIR/codon_full/$stem.codon.full.fasta"
  codon="$OUTPUT_DIR/codon_aligned/$stem.codon.fasta"
  log="$OUTPUT_DIR/logs/$stem.alignment.log"

  if [[ ! -s "$cds" ]]; then
    echo "Missing CDS FASTA: $cds" >&2
    exit 1
  fi
  if [[ "$RESUME" == "1" && -s "$aligned" && -s "$trimmed" && -s "$codon_full" && -s "$codon" ]]; then
    printf '%s\tSKIP_COMPLETE\n' "$stem" >> "$OUTPUT_DIR/alignment_status.tsv"
    skipped=$((skipped + 1))
    continue
  fi

  mafft --auto --thread "$THREADS" "$protein" > "$aligned" 2> "$log"
  # PAL2NAL must receive the complete protein alignment. A protein alignment
  # with columns already removed by trimAl is no longer consistent with the
  # complete CDS and causes a false peptide/nucleotide inconsistency error.
  pal2nal.pl "$aligned" "$cds" -codontable 11 -output fasta > "$codon_full" 2>> "$log"
  trimal -in "$aligned" -out "$trimmed" -automated1 -keepheader >> "$log" 2>&1
  # Apply exactly the same protein-column mask while back-translating the CDS,
  # so every retained amino-acid column corresponds to one retained codon.
  trimal -in "$aligned" -out "$codon" -automated1 -keepheader \
    -backtrans "$cds" -ignorestopcodon >> "$log" 2>&1

  if [[ ! -s "$codon_full" || ! -s "$trimmed" || ! -s "$codon" ]]; then
    echo "Empty output for family: $stem" >&2
    exit 1
  fi
  printf '%s\tPASS\n' "$stem" >> "$OUTPUT_DIR/alignment_status.tsv"
  processed=$((processed + 1))
done

protein_count=$(find "$OUTPUT_DIR/protein_trimmed" -maxdepth 1 -name '*.trim.faa' -size +0c | wc -l)
codon_count=$(find "$OUTPUT_DIR/codon_aligned" -maxdepth 1 -name '*.codon.fasta' -size +0c | wc -l)
codon_full_count=$(find "$OUTPUT_DIR/codon_full" -maxdepth 1 -name '*.codon.full.fasta' -size +0c | wc -l)
echo "Newly processed families: $processed"
echo "Previously complete families skipped: $skipped"
echo "Trimmed protein alignments present: $protein_count"
echo "Codon alignments present: $codon_count"
echo "Full PAL2NAL codon alignments present: $codon_full_count"
[[ "$protein_count" -eq 112 && "$codon_count" -eq 112 && "$codon_full_count" -eq 112 ]]
