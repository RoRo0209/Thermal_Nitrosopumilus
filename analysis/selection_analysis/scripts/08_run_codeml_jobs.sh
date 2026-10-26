#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 PAML_JOB_DIR [PARALLEL_JOBS]" >&2
  exit 2
fi

ROOT="$(realpath "$1")"
JOBS="${2:-8}"

if ! command -v codeml >/dev/null 2>&1; then
  echo "codeml not found in PATH" >&2
  exit 1
fi

find "$ROOT" -name codeml.ctl -print0 | \
  xargs -0 -r -n1 -P "$JOBS" bash -c '
    ctl="$1"
    dir="$(dirname "$ctl")"
    cd "$dir"
    if codeml codeml.ctl > codeml.stdout 2> codeml.stderr; then
      printf "PASS\t%s\n" "$dir"
    else
      printf "FAIL\t%s\n" "$dir" >&2
      exit 1
    fi
  ' _

echo "All codeml jobs completed under: $ROOT"
