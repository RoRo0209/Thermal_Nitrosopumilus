#!/usr/bin/env python3
"""Prepare a targeted 30 x 2 multi-start codeml grid for iscU2 only."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


OMEGA_STARTS = (0.005, 0.02, 0.05, 0.2, 1.0, 5.0)
KAPPA_STARTS = (0.25, 0.75, 2.0, 5.0, 12.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-family", type=Path, required=True)
    parser.add_argument("--original-summary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def tag(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def ctl_text(treefile: str, omega: float, kappa: float) -> str:
    return f"""seqfile = ../alignment.phy
treefile = ../{treefile}
outfile = results.txt

noisy = 3
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = 2
NSsites = 0
icode = 0
Mgene = 0
fix_kappa = 0
kappa = {kappa:g}
fix_omega = 0
omega = {omega:g}
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 10
getSE = 0
RateAncestor = 1
Small_Diff = 1e-6
cleandata = 1
method = 0
"""


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    jobs = args.outdir / "jobs" / "iscU2"
    jobs.mkdir(parents=True, exist_ok=True)

    required = {
        "alignment.phy": args.source_family / "alignment.phy",
        "null_H_equals_L.tree": args.source_family / "null_H_equals_L.tree",
        "alternative_H_L_O.tree": args.source_family / "alternative_H_L_O.tree",
    }
    for target_name, source in required.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, jobs / target_name)

    if not args.original_summary.is_file():
        raise FileNotFoundError(args.original_summary)
    inputs = args.outdir / "inputs"
    inputs.mkdir(exist_ok=True)
    shutil.copy2(
        args.original_summary,
        inputs / "direct_HighT_LowT_omega_LRT_original.tsv",
    )

    rows = []
    for model, treefile in (
        ("null", "null_H_equals_L.tree"),
        ("alternative", "alternative_H_L_O.tree"),
    ):
        for omega in OMEGA_STARTS:
            for kappa in KAPPA_STARTS:
                run_id = f"{model}_o{tag(omega)}_k{tag(kappa)}"
                run_dir = jobs / run_id
                run_dir.mkdir(exist_ok=True)
                (run_dir / "codeml.ctl").write_text(
                    ctl_text(treefile, omega, kappa),
                    encoding="ascii",
                    newline="\n",
                )
                rows.append(
                    {
                        "family": "iscU2",
                        "model": model,
                        "run_id": run_id,
                        "start_omega": omega,
                        "start_kappa": kappa,
                        "treefile": treefile,
                        "ctl": str(run_dir / "codeml.ctl"),
                    }
                )

    with (args.outdir / "iscU2_multistart_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"iscU2 targeted runs prepared: {len(rows)}")
    print(f"Null runs: {sum(row['model'] == 'null' for row in rows)}")
    print(f"Alternative runs: {sum(row['model'] == 'alternative' for row in rows)}")
    print("Omega starts:", ",".join(map(str, OMEGA_STARTS)))
    print("Kappa starts:", ",".join(map(str, KAPPA_STARTS)))


if __name__ == "__main__":
    main()
