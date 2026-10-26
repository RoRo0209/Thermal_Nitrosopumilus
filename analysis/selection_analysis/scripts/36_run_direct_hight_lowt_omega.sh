#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 DIRECT_JOB_DIR [PARALLEL_JOBS]" >&2
  exit 2
fi

ROOT="$(realpath "$1")"
JOBS="${2:-16}"

if ! command -v codeml >/dev/null 2>&1; then
  echo "Required executable not found: codeml" >&2
  exit 1
fi

find "$ROOT" -name codeml.ctl -print0 | \
  xargs -0 -r -n1 -P "$JOBS" bash -c '
    ctl="$1"
    dir="$(dirname "$ctl")"
    cd "$dir"
    if [[ -s results.txt ]] && grep -q "Time used:" results.txt; then
      printf "SKIP_COMPLETE\t%s\n" "$dir"
    elif codeml codeml.ctl > codeml.stdout 2> codeml.stderr; then
      if grep -q "Time used:" results.txt; then
        printf "PASS\t%s\n" "$dir"
      else
        printf "FAIL_NO_COMPLETION_MARKER\t%s\n" "$dir" >&2
        exit 1
      fi
    else
      printf "FAIL\t%s\n" "$dir" >&2
      exit 1
    fi
  ' _

echo "All direct HighT-LowT codeml jobs completed under: $ROOT"
