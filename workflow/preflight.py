#!/usr/bin/env python3
"""Check released file integrity, or list input gaps without running analysis."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def code_check():
    failures = []
    count = 0
    for line in (ROOT / "CHECKSUMS.sha256").read_text(encoding="utf-8").splitlines():
        digest, rel = line.split("  ", 1)
        path = ROOT / rel
        count += 1
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            failures.append("Missing or changed: " + rel)
    for path in ROOT.rglob("*.py"):
        try:
            ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        except SyntaxError as exc:
            failures.append(str(exc))
    for message in failures:
        print(message)
    print(f"Checked {count} release files; Python syntax checked; failures={len(failures)}.")
    print("This check does not execute scientific analyses or verify manuscript results.")
    return 1 if failures else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--code-only", action="store_true")
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--module", help="Module directory name, e.g. ALE_2 or 60weeks")
    args = parser.parse_args()
    if args.code_only:
        return code_check()
    if args.workspace is None:
        parser.error("Choose --code-only or --workspace PATH [--module NAME]")
    stages = json.loads((ROOT / "workflow/stages.json").read_text(encoding="utf-8"))
    selected = [s for s in stages if not args.module or s["module"] == args.module]
    if not selected:
        parser.error("Unknown module")
    print("Input inventory only; free-text/CLI requirements need manual reconciliation.")
    for stage in selected:
        print(f'\n{stage["id"]}: {stage["status"]}')
        for spec in stage["inputs"]:
            if spec.startswith("env:"):
                key = spec[4:]
                value = os.environ.get(key)
                print(f"  {key}: " + ("SET (schema/content not verified)" if value else "MISSING"))
            elif any(ch in spec for ch in (" ", "<", ">", ";")) or spec.startswith("arg") or "/" not in spec:
                print("  MANUAL: " + spec)
            else:
                path = args.workspace / spec
                status = "PRESENT (not validated)" if path.exists() else "MISSING or generated later"
                print(f"  {status}: {spec}")
    print("\nINCOMPLETE: authorized data, exact provenance, versions and manual stage requirements remain necessary.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
