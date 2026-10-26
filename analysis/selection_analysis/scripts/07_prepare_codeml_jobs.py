#!/usr/bin/env python3
"""Validate codon alignments and create replicated codeml foreground jobs."""

from __future__ import annotations

import argparse
import csv
import shutil
from collections import Counter
from pathlib import Path


STOPS = {"TAA", "TAG", "TGA"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--codon-dir", type=Path, required=True)
    p.add_argument("--tree-dir", type=Path, required=True)
    p.add_argument("--family-qc", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--min-codons", type=int, default=50)
    p.add_argument("--max-gap-fraction", type=float, default=0.50)
    return p.parse_args()


def read_fasta(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    result[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        result[name] = "".join(chunks).upper()
    return result


def write_phylip(path: Path, sequences: dict[str, str]) -> None:
    length = len(next(iter(sequences.values())))
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(f" {len(sequences)} {length}\n")
        for name, sequence in sequences.items():
            handle.write(f"{name}  {sequence}\n")


def ctl_text(treefile: str, model: int, nssites: int, fix_omega: int, omega: float) -> str:
    return f"""seqfile = ../../alignment.phy
treefile = ../{treefile}
outfile = results.txt

noisy = 3
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = {model}
NSsites = {nssites}
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = {fix_omega}
omega = {omega}
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 10
getSE = 0
RateAncestor = 0
Small_Diff = 1e-6
cleandata = 1
method = 0
"""


def make_model(ecotype_dir: Path, name: str, treefile: str, model: int, nssites: int, fix_omega: int, starts: list[float]):
    for i, omega in enumerate(starts, start=1):
        run_dir = ecotype_dir / f"{name}_rep{i}"
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "codeml.ctl").write_text(
            ctl_text(treefile, model, nssites, fix_omega, omega), encoding="ascii", newline="\n"
        )


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    family_map = {}
    with args.family_qc.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            family_map[row["file_stem"]] = row["split_family"]

    rows = []
    for codon_path in sorted(args.codon_dir.glob("*.codon.fasta")):
        stem = codon_path.name.removesuffix(".codon.fasta")
        sequences = read_fasta(codon_path)
        lengths = Counter(map(len, sequences.values()))
        equal_length = len(lengths) == 1
        alignment_length = next(iter(lengths)) if equal_length and lengths else 0
        frame_ok = alignment_length % 3 == 0
        codons = alignment_length // 3 if frame_ok else 0
        gap_fraction = (
            sum(seq.count("-") for seq in sequences.values()) /
            (len(sequences) * alignment_length)
            if sequences and alignment_length else 1.0
        )
        internal_stop_sequences = 0
        for sequence in sequences.values():
            ungapped_codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
            if any(codon in STOPS for codon in ungapped_codons[:-1]):
                internal_stop_sequences += 1
        ids = set(sequences)

        tree_status = {}
        for ecotype in ("HighT", "LowT"):
            for scope in ("stem", "clade"):
                tree_path = args.tree_dir / f"foreground_{ecotype}_{scope}" / f"{stem}.nwk"
                tree_status[f"{ecotype}_{scope}"] = tree_path.is_file() and tree_path.stat().st_size > 0

        eligible = (
            len(sequences) >= 4 and equal_length and frame_ok and codons >= args.min_codons
            and gap_fraction <= args.max_gap_fraction and internal_stop_sequences == 0
            and all(tree_status.values()) and stem in family_map
        )
        rows.append({
            "split_family": family_map.get(stem, ""),
            "file_stem": stem,
            "n_sequences": len(sequences),
            "alignment_nt": alignment_length,
            "alignment_codons": codons,
            "equal_length": int(equal_length),
            "frame_ok": int(frame_ok),
            "gap_fraction": f"{gap_fraction:.6g}",
            "internal_stop_sequences": internal_stop_sequences,
            "HighT_stem_tree_found": int(tree_status["HighT_stem"]),
            "LowT_stem_tree_found": int(tree_status["LowT_stem"]),
            "HighT_clade_tree_found": int(tree_status["HighT_clade"]),
            "LowT_clade_tree_found": int(tree_status["LowT_clade"]),
            "codeml_eligible": int(eligible),
        })
        if not eligible:
            continue

        family_dir = args.outdir / stem
        family_dir.mkdir(exist_ok=True)
        write_phylip(family_dir / "alignment.phy", sequences)
        for ecotype in ("HighT", "LowT"):
            ecotype_dir = family_dir / ecotype
            ecotype_dir.mkdir(exist_ok=True)
            # PAML 4.10 expects the PHYLIP tree header "n_taxa n_trees".
            # Keep generic .nwk inputs as plain Newick and add the header only
            # to the codeml-specific copies.
            for source, destination in (
                (
                    args.tree_dir / f"foreground_{ecotype}_clade" / f"{stem}.nwk",
                    ecotype_dir / "branch_foreground.tree",
                ),
                (
                    args.tree_dir / f"foreground_{ecotype}_stem" / f"{stem}.nwk",
                    ecotype_dir / "branchsite_foreground.tree",
                ),
            ):
                newick = source.read_text(encoding="utf-8").strip()
                destination.write_text(
                    f"{len(sequences)} 1\n{newick}\n", encoding="ascii", newline="\n"
                )
            make_model(ecotype_dir, "branch_null", "branch_foreground.tree", 0, 0, 0, [0.2])
            make_model(ecotype_dir, "branch_alt", "branch_foreground.tree", 2, 0, 0, [0.2, 1.0, 2.0])
            make_model(ecotype_dir, "branchsite_null", "branchsite_foreground.tree", 2, 2, 1, [1.0])
            make_model(ecotype_dir, "branchsite_alt", "branchsite_foreground.tree", 2, 2, 0, [0.5, 1.5, 3.0])

    fields = list(rows[0]) if rows else []
    with (args.outdir / "codon_alignment_qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Codon alignments inspected: {len(rows)}")
    print(f"codeml-ready families: {sum(int(row['codeml_eligible']) for row in rows)}")


if __name__ == "__main__":
    main()
