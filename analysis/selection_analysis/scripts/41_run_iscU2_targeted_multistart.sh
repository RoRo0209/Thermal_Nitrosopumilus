#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:?Usage: $0 TARGET_ROOT [PARALLEL_CODEML]}"
PARALLEL_CODEML="${2:-${PARALLEL_CODEML:-16}}"
JOBS="$ROOT/jobs/iscU2"

if ! command -v codeml >/dev/null 2>&1; then
  echo "Required executable not found: codeml" >&2
  exit 1
fi

find "$JOBS" -mindepth 2 -maxdepth 2 -name codeml.ctl -print0 |
  xargs -0 -r -n1 -P "$PARALLEL_CODEML" bash -c '
    ctl="$1"
    run_dir=$(dirname "$ctl")
    if [[ -s "$run_dir/results.txt" ]] &&
       grep -q "Time used:" "$run_dir/results.txt"; then
      printf "SKIP\t%s\n" "$run_dir"
      exit 0
    fi
    (
      cd "$run_dir"
      codeml codeml.ctl > codeml.stdout 2> codeml.stderr
    )
    if [[ -s "$run_dir/results.txt" ]] &&
       grep -q "Time used:" "$run_dir/results.txt"; then
      printf "PASS\t%s\n" "$run_dir"
    else
      printf "FAIL\t%s\n" "$run_dir"
      exit 1
    fi
  ' _

echo "All iscU2 targeted multi-start jobs completed under: $JOBS"
