#!/usr/bin/env python
"""
Summarize HighT-vs-LowT allele-frequency differentiation in the SPOT01 C2/amoABC region.

Usage:
  python amoabc_site_differentiation.py \
    --vcf amoabc_workstation_out/amoabc.norm.vcf.gz \
    --groups amoabc_workstation_out/sample_groups.tsv \
    --ref refs/LowT_representative.fna \
    --gff refs/LowT_representative.gff \
    --outdir amoabc_site_diff

The script uses FORMAT/AD when available, so it can summarize allele fractions
from metagenomic read recruitments instead of relying only on hard genotypes.
"""

from __future__ import annotations

# Release path adapter: data and outputs live in a separate user workspace.
import os as _release_os
from pathlib import Path as _ReleasePath
if not _release_os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the analysis workspace")
_release_root = _ReleasePath(_release_os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
_release_os.chdir(_release_root / 'tajima')


import argparse
import gzip
import math
import os
from pathlib import Path


# Reference-coordinate loci must be supplied with the matched reference.
import csv
with open(_release_os.environ["THERMAL_AMO_LOCI"], encoding="utf-8") as _loci_handle:
    _loci = list(csv.DictReader(_loci_handle, delimiter="\t"))
AMO_GENES = {row["gene"]: (row["gene_id"], int(row["start"]), int(row["end"]), row["strand"]) for row in _loci}
if set(AMO_GENES) != {"amoA", "amoB", "amoC"} or len(_loci) != 3:
    raise ValueError("THERMAL_AMO_LOCI must define amoA, amoB and amoC exactly once")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def read_groups(path: Path) -> dict[str, str]:
    groups = {}
    with open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] == "sample":
                continue
            groups[parts[0]] = parts[1]
    return groups


def resolve_clade(sample: str, groups: dict[str, str]) -> str | None:
    if sample in groups:
        return groups[sample]
    base = os.path.basename(sample)
    base = base.replace(".bam", "").replace(".cram", "")
    for suffix in ("_HighT", "_LowT", ".HighT", ".LowT"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    if base in groups:
        return groups[base]
    matches = [clade for sample_id, clade in groups.items() if sample_id in sample or sample_id in base]
    if len(set(matches)) == 1:
        return matches[0]
    return None


def locate_amo_gene(pos: int):
    for gene, (_, start, end, strand) in AMO_GENES.items():
        if start <= pos <= end:
            return gene, start, end, strand
    return "outside_amoABC", None, None, None


def classify_variant(chrom_seq: str, pos: int, ref: str, alt: str):
    gene, start, end, strand = locate_amo_gene(pos)
    if gene == "outside_amoABC" or len(ref) != 1 or len(alt) != 1:
        return gene, "", "", "", "", ""

    assert start is not None and end is not None and strand is not None
    if strand == "+":
        offset = pos - start
        codon_start = start + (offset // 3) * 3
        codon = chrom_seq[codon_start - 1: codon_start + 2]
        idx = offset % 3
        alt_codon = codon[:idx] + alt.upper() + codon[idx + 1:]
        codon_pos = offset // 3 + 1
    else:
        offset = end - pos
        genomic_codon_end = end - (offset // 3) * 3
        genomic_codon_start = genomic_codon_end - 2
        genomic_codon = chrom_seq[genomic_codon_start - 1: genomic_codon_end]
        codon = revcomp(genomic_codon)
        idx = offset % 3
        coding_alt = revcomp(alt.upper())
        alt_codon = codon[:idx] + coding_alt + codon[idx + 1:]
        codon_pos = offset // 3 + 1

    aa_ref = CODON_TABLE.get(codon.upper(), "X")
    aa_alt = CODON_TABLE.get(alt_codon.upper(), "X")
    effect = "synonymous" if aa_ref == aa_alt else "nonsynonymous"
    return gene, str(codon_pos), codon.upper(), alt_codon.upper(), aa_ref, aa_alt, effect


def parse_ad(sample_field: str, fmt_keys: list[str], alt_index: int = 1):
    values = sample_field.split(":")
    data = dict(zip(fmt_keys, values))
    ad = data.get("AD")
    dp = data.get("DP")
    if ad and ad not in (".", ""):
        parts = [int(x) if x not in (".", "") else 0 for x in ad.split(",")]
        ref_depth = parts[0] if len(parts) > 0 else 0
        alt_depth = parts[alt_index] if len(parts) > alt_index else 0
        total = sum(parts)
        return ref_depth, alt_depth, total
    if dp and dp not in (".", ""):
        return 0, 0, int(dp)
    return 0, 0, 0


def hudson_fst_from_counts(ref_h, alt_h, ref_l, alt_l):
    n1 = ref_h + alt_h
    n2 = ref_l + alt_l
    if n1 <= 1 or n2 <= 1:
        return math.nan
    p1 = alt_h / n1
    p2 = alt_l / n2
    h1 = 2 * p1 * (1 - p1) * n1 / (n1 - 1)
    h2 = 2 * p2 * (1 - p2) * n2 / (n2 - 1)
    between = (p1 - p2) ** 2
    denom = between + h1 / 2 + h2 / 2
    if denom <= 0:
        return math.nan
    return between / denom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--groups", required=True, type=Path)
    ap.add_argument("--ref", required=True, type=Path)
    ap.add_argument("--gff", type=Path, help="Kept for provenance; amoABC coordinates are embedded.")
    ap.add_argument("--outdir", default=Path("amoabc_site_diff"), type=Path)
    ap.add_argument("--min-depth-per-pop", default=20, type=int)
    ap.add_argument("--fixed-threshold", default=0.95, type=float)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    groups = read_groups(args.groups)
    fasta = read_fasta(args.ref)

    records = []
    samples = []
    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos_s, vid, ref, alt_s, qual, filt, info, fmt = parts[:9]
            pos = int(pos_s)
            alts = alt_s.split(",")
            # bcftools norm -m-any should leave one ALT; keep first if not.
            alt = alts[0]
            fmt_keys = fmt.split(":")
            chrom_seq = fasta.get(chrom)
            if chrom_seq is None:
                continue

            counts = {
                "HighT": {"ref": 0, "alt": 0, "depth": 0, "n": 0},
                "LowT": {"ref": 0, "alt": 0, "depth": 0, "n": 0},
            }
            for sample, sample_field in zip(samples, parts[9:]):
                clade = resolve_clade(sample, groups)
                if clade not in counts:
                    continue
                rd, ad, td = parse_ad(sample_field, fmt_keys)
                counts[clade]["ref"] += rd
                counts[clade]["alt"] += ad
                counts[clade]["depth"] += td
                if td > 0:
                    counts[clade]["n"] += 1

            h = counts["HighT"]
            l = counts["LowT"]
            if h["depth"] < args.min_depth_per_pop or l["depth"] < args.min_depth_per_pop:
                continue
            h_af = h["alt"] / h["depth"] if h["depth"] else math.nan
            l_af = l["alt"] / l["depth"] if l["depth"] else math.nan
            delta = h_af - l_af
            abs_delta = abs(delta)
            fixed_like = (
                (h_af >= args.fixed_threshold and l_af <= 1 - args.fixed_threshold)
                or (l_af >= args.fixed_threshold and h_af <= 1 - args.fixed_threshold)
            )
            gene_info = classify_variant(chrom_seq, pos, ref, alt)
            gene = gene_info[0]
            if gene == "outside_amoABC":
                codon_pos = ref_codon = alt_codon = aa_ref = aa_alt = effect = ""
            else:
                _, codon_pos, ref_codon, alt_codon, aa_ref, aa_alt, effect = gene_info

            records.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": gene,
                "codon_pos": codon_pos,
                "ref_codon": ref_codon,
                "alt_codon": alt_codon,
                "aa_ref": aa_ref,
                "aa_alt": aa_alt,
                "effect": effect,
                "HighT_ref_depth": h["ref"],
                "HighT_alt_depth": h["alt"],
                "HighT_total_depth": h["depth"],
                "HighT_n_samples": h["n"],
                "HighT_alt_freq": h_af,
                "LowT_ref_depth": l["ref"],
                "LowT_alt_depth": l["alt"],
                "LowT_total_depth": l["depth"],
                "LowT_n_samples": l["n"],
                "LowT_alt_freq": l_af,
                "delta_alt_freq_HighT_minus_LowT": delta,
                "abs_delta_alt_freq": abs_delta,
                "fixed_like": fixed_like,
                "site_fst_proxy": hudson_fst_from_counts(h["ref"], h["alt"], l["ref"], l["alt"]),
            })

    out_tsv = args.outdir / "amoabc_site_differentiation.tsv"
    fields = [
        "chrom", "pos", "ref", "alt", "gene", "codon_pos", "ref_codon", "alt_codon",
        "aa_ref", "aa_alt", "effect", "HighT_ref_depth", "HighT_alt_depth",
        "HighT_total_depth", "HighT_n_samples", "HighT_alt_freq", "LowT_ref_depth",
        "LowT_alt_depth", "LowT_total_depth", "LowT_n_samples", "LowT_alt_freq",
        "delta_alt_freq_HighT_minus_LowT", "abs_delta_alt_freq", "fixed_like",
        "site_fst_proxy",
    ]
    with open(out_tsv, "wt") as out:
        out.write("\t".join(fields) + "\n")
        for r in sorted(records, key=lambda x: (x["chrom"], x["pos"])):
            out.write("\t".join(str(r[f]) for f in fields) + "\n")

    amo_records = [r for r in records if r["gene"] in {"amoA", "amoB", "amoC"}]
    nonsyn = [r for r in amo_records if r["effect"] == "nonsynonymous"]
    fixed = [r for r in amo_records if r["fixed_like"]]
    high_delta = [r for r in amo_records if r["abs_delta_alt_freq"] >= 0.5]

    with open(args.outdir / "amoabc_site_differentiation_summary.txt", "wt") as out:
        out.write("amoABC site differentiation summary\n")
        out.write("===================================\n")
        out.write(f"VCF samples parsed: {len(samples)}\n")
        out.write(f"Sites passing depth filter in C2 region: {len(records)}\n")
        out.write(f"Sites inside amoABC: {len(amo_records)}\n")
        out.write(f"Nonsynonymous sites inside amoABC: {len(nonsyn)}\n")
        out.write(f"Fixed-like sites inside amoABC: {len(fixed)}\n")
        out.write(f"Sites with abs(delta AF) >= 0.5 inside amoABC: {len(high_delta)}\n")
        out.write("\nTop amoABC sites by abs(delta AF):\n")
        for r in sorted(amo_records, key=lambda x: x["abs_delta_alt_freq"], reverse=True)[:20]:
            out.write(
                f"{r['gene']}\t{r['pos']}\t{r['ref']}>{r['alt']}\t"
                f"HighT_AF={r['HighT_alt_freq']:.3f}\tLowT_AF={r['LowT_alt_freq']:.3f}\t"
                f"delta={r['delta_alt_freq_HighT_minus_LowT']:.3f}\t"
                f"{r['effect']}\t{r['aa_ref']}{r['codon_pos']}{r['aa_alt']}\n"
            )

    print(f"Wrote {out_tsv}")
    print(f"Wrote {args.outdir / 'amoabc_site_differentiation_summary.txt'}")


if __name__ == "__main__":
    main()
