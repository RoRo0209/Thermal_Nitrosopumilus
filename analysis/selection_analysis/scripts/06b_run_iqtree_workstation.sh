#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 PROTEIN_TRIMMED_DIR GENE_TREE_OUTPUT_DIR" >&2
  exit 2
fi

PROTEIN_DIR="$(realpath "$1")"
OUTPUT_DIR="$2"
FAMILY_JOBS="${FAMILY_JOBS:-4}"
THREADS_PER_JOB="${THREADS_PER_JOB:-4}"
IQTREE_BOOTSTRAPS="${IQTREE_BOOTSTRAPS:-0}"
mkdir -p "$OUTPUT_DIR/logs"

IQTREE=""
for candidate in iqtree3 iqtree2 iqtree; do
  if command -v "$candidate" >/dev/null 2>&1; then
    IQTREE="$candidate"
    break
  fi
done
if [[ -z "$IQTREE" ]]; then
  echo "IQ-TREE executable not found (tried iqtree3, iqtree2, iqtree)" >&2
  exit 1
fi

export OUTPUT_DIR THREADS_PER_JOB IQTREE_BOOTSTRAPS IQTREE
find "$PROTEIN_DIR" -maxdepth 1 -type f -name '*.trim.faa' -print0 | \
  xargs -0 -r -n1 -P "$FAMILY_JOBS" bash -c '
    set -euo pipefail
    aln="$1"
    stem="$(basename "$aln" .trim.faa)"
    prefix="$OUTPUT_DIR/$stem"
    log="$OUTPUT_DIR/logs/$stem.iqtree.log"
    if [[ -s "$prefix.treefile" ]]; then
      echo "SKIP_COMPLETE $stem"
      exit 0
    fi
    args=(-s "$aln" -m MFP -T "$THREADS_PER_JOB" --prefix "$prefix")
    if [[ "$IQTREE_BOOTSTRAPS" -gt 0 ]]; then
      args+=(-B "$IQTREE_BOOTSTRAPS")
    fi
    "$IQTREE" "${args[@]}" > "$log" 2>&1
    echo "PASS $stem"
  ' _

expected=$(find "$PROTEIN_DIR" -maxdepth 1 -name '*.trim.faa' -size +0c | wc -l)
complete=$(find "$OUTPUT_DIR" -maxdepth 1 -name '*.treefile' -size +0c | wc -l)
echo "IQ-TREE executable: $IQTREE"
echo "Expected/completed gene trees: $expected/$complete"
[[ "$complete" -eq "$expected" ]]

