#!/usr/bin/env python3
"""Extract branch events from every family in ONE original ALEml_undated run.

The species-tree S: record must equal the record verified by
verify_ale_mapping.R. No other run, precomputed mapping, or R index is used.
"""
import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path


def parse_reconciliation(path):
    species, records = None, {}
    with path.open(encoding="utf-8-sig") as handle:
        header = handle.readline()
        if not header.startswith("#ALEml_undated"):
            raise ValueError(f"Not an ALEml_undated output: {path.name}")
        for line in handle:
            if line.startswith("S:"):
                if species is not None:
                    raise ValueError(f"Repeated S: record: {path.name}")
                species = "".join(line[2:].split())
            if not line.startswith(("S_internal_branch", "S_terminal_branch")):
                continue
            fields = line.split()
            if len(fields) < 6 or fields[1] in records:
                raise ValueError(f"Malformed/duplicate branch record: {path.name}")
            branch = fields[1]
            if not re.fullmatch(r"\d+|[^\s()]+\(\d+\)", branch):
                raise ValueError(f"Unrecognized ALE branch: {branch}")
            if (fields[0] == "S_terminal_branch") != ("(" in branch):
                raise ValueError(f"Branch kind/ID mismatch: {path.name}")
            try:
                values = tuple(Decimal(x) for x in fields[2:6])
            except InvalidOperation as exc:
                raise ValueError(f"Non-numeric event: {path.name}") from exc
            if any(not x.is_finite() or x < 0 for x in values):
                raise ValueError(f"Invalid event value: {path.name}")
            records[branch] = values
    if species is None or not records:
        raise ValueError(f"Missing S: tree or branch events: {path.name}")
    return species, records


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fields)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uml-dir", required=True, type=Path)
    parser.add_argument("--verified-dir", required=True, type=Path,
                        help="Directory written by verify_ale_mapping.R")
    parser.add_argument("--expected-families", required=True, type=int,
                        help="Expected family count from this run's provenance")
    args = parser.parse_args()
    files = sorted(args.uml_dir.rglob("*.uml_rec"))
    if args.expected_families <= 0 or len(files) != args.expected_families:
        parser.error(f"Expected {args.expected_families} original family outputs; found {len(files)}")
    out = args.verified_dir
    reference = "".join((out / "ale_numbered_species_tree.nwk").read_text(encoding="utf-8").split())
    with (out / "crowns.tsv").open(encoding="utf-8", newline="") as handle:
        crowns = list(csv.DictReader(handle, delimiter="\t"))
    if len(crowns) != 2 or {x["group"] for x in crowns} != {"HighT", "LowT"}:
        raise ValueError("Expected exactly one verified crown per group")
    branches = {x["ALE_ID"] for x in crowns}
    if len(branches) != 2 or any(not x.isdigit() for x in branches):
        raise ValueError("Expected distinct internal ALE crown IDs")
    sums = defaultdict(lambda: [Decimal(0)] * 4)
    selected = {b: [] for b in branches}
    with (out / "node_mapping.tsv").open(encoding="utf-8", newline="") as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    expected_records = {
        f'{r["descendant_ids"]}({r["ALE_ID"]})' if r["isTip"] == "TRUE" else r["ALE_ID"]
        for r in mapping
    }
    if not mapping or len(expected_records) != len(mapping) or len({r["ALE_ID"] for r in mapping}) != len(mapping):
        raise ValueError("Invalid or duplicate verified node mapping")
    for crown in crowns:
        matches = [r for r in mapping if r["ALE_ID"] == crown["ALE_ID"] and r["isTip"] == "FALSE"]
        if len(matches) != 1 or matches[0]["descendant_ids"] != crown["descendant_ids"]:
            raise ValueError("Crown membership differs from verified node mapping")
    seen, sources = set(), []
    for path in files:
        species, records = parse_reconciliation(path)
        if species != reference:
            raise ValueError(f"S: tree/numbering differs from verified reference: {path.name}")
        if not branches <= records.keys():
            raise ValueError(f"Crown missing from branch records: {path.name}")
        if set(records) != expected_records:
            raise ValueError(f"Branch coverage differs: {path.name}")
        suffix = ".ale.uml_rec"
        if not path.name.endswith(suffix):
            raise ValueError(f"Unrecognized family output suffix: {path.name}")
        gene = path.name[:-len(suffix)]
        if gene in seen:
            raise ValueError(f"Duplicate family output: {gene}")
        seen.add(gene)
        for branch, events in records.items():
            sums[branch] = [a + b for a, b in zip(sums[branch], events)]
        for branch in branches:
            selected[branch].append((gene, *records[branch]))
        sources.append({"file": str(path.relative_to(args.uml_dir)),
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    write_tsv(out / "DTLO_table_rgy.txt", ["Branch", "D", "T", "L", "O"],
              ((b, *sums[b]) for b in sorted(sums)))
    for branch in sorted(branches):
        write_tsv(out / f"branch{branch}_gene_DTLO.tsv", ["gene", "D", "T", "L", "O"], selected[branch])
    output_names = ["crowns.tsv", "node_mapping.tsv", "ale_numbered_species_tree.nwk", "DTLO_table_rgy.txt"]
    output_names += [f"branch{branch}_gene_DTLO.tsv" for branch in sorted(branches)]
    receipt = {"n_families": len(seen), "n_branches": len(sums), "crowns": crowns,
               "outputs_md5": {name: hashlib.md5((out / name).read_bytes()).hexdigest() for name in output_names},
               "species_tree_sha256": hashlib.sha256(reference.encode()).hexdigest(),
               "sources": sources,
               "event_definition": "Expected counts from columns D/T/L/O; all families retained; no threshold."}
    (out / "same_run_extraction.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(f"Extracted {len(seen)} families and {len(sums)} branches from one matching numbered species tree.")


if __name__ == "__main__":
    main()
