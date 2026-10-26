#!/usr/bin/env python3
"""Extract exact-run Prokka CDS/proteins for Roary soft-core families and QC them."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}
TABLE11_STARTS = {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--members", type=Path, required=True)
    p.add_argument("--ffn-dir", type=Path, required=True)
    p.add_argument("--faa-dir", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--min-hight", type=int, default=10)
    p.add_argument("--min-lowt", type=int, default=20)
    return p.parse_args()


def read_fasta(path: Path):
    header = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks).upper()


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def normalize_locus(header: str, kind: str) -> str:
    token = header.split()[0]
    if kind == "ffn":
        return token.split("$", 1)[-1]
    return token.split("__", 1)[-1]


def load_sequences(directory: Path, suffix: str, kind: str):
    sequences: dict[str, str] = {}
    source: dict[str, str] = {}
    duplicates: set[str] = set()
    for path in sorted(directory.glob(f"*.{suffix}")):
        for header, seq in read_fasta(path):
            locus = normalize_locus(header, kind)
            if locus in sequences:
                duplicates.add(locus)
            sequences[locus] = seq
            source[locus] = path.stem
    if duplicates:
        examples = ", ".join(sorted(duplicates)[:10])
        raise ValueError(f"Duplicate {kind} locus tags ({len(duplicates)}): {examples}")
    return sequences, source


def translate_cds(cds: str) -> str:
    aa = []
    usable = len(cds) - (len(cds) % 3)
    for i in range(0, usable, 3):
        aa.append(CODON_TABLE.get(cds[i:i + 3], "X"))
    return "".join(aa)


def translation_status(cds: str, protein: str):
    translated = translate_cds(cds)
    translated_no_stop = translated[:-1] if translated.endswith("*") else translated
    protein_no_stop = protein[:-1] if protein.endswith("*") else protein
    internal_stops = translated_no_stop.count("*")
    if translated_no_stop == protein_no_stop:
        status = "exact"
    elif (
        len(translated_no_stop) == len(protein_no_stop)
        and translated_no_stop[1:] == protein_no_stop[1:]
        and protein_no_stop.startswith("M")
        and cds[:3] in TABLE11_STARTS
    ):
        status = "exact_table11_start"
    else:
        status = "mismatch"
    mismatch_positions = [
        i + 1 for i, (a, b) in enumerate(zip(translated_no_stop, protein_no_stop)) if a != b
    ]
    length_delta = len(translated_no_stop) - len(protein_no_stop)
    mismatch_count = len(mismatch_positions) + abs(length_delta)
    first_mismatch = mismatch_positions[0] if mismatch_positions else (min(len(translated_no_stop), len(protein_no_stop)) + 1 if length_delta else "")
    first_pair = ""
    if mismatch_positions:
        i = mismatch_positions[0] - 1
        first_pair = f"{translated_no_stop[i]}>{protein_no_stop[i]}"
    return translated, status, internal_stops, mismatch_count, first_mismatch, first_pair, length_delta


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def tree_id_from_locus(locus: str) -> str:
    match = re.match(r"^(.+)_\d+$", locus)
    if not match:
        raise ValueError(f"Cannot derive Prokka/tree prefix from locus tag: {locus}")
    return match.group(1)


def write_tsv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cds_dir = args.outdir / "raw_cds"
    protein_dir = args.outdir / "raw_proteins"
    clean_cds_dir = args.outdir / "qc_cds"
    clean_protein_dir = args.outdir / "qc_proteins"
    cds_dir.mkdir(exist_ok=True)
    protein_dir.mkdir(exist_ok=True)
    clean_cds_dir.mkdir(exist_ok=True)
    clean_protein_dir.mkdir(exist_ok=True)

    families: dict[str, list[dict[str, str]]] = defaultdict(list)
    with args.members.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            families[row["split_family"]].append(row)

    ffn, ffn_source = load_sequences(args.ffn_dir, "ffn", "ffn")
    faa, faa_source = load_sequences(args.faa_dir, "faa", "faa")

    sequence_rows: list[dict] = []
    family_rows: list[dict] = []
    eligible_families: list[str] = []
    genome_to_tree_id: dict[str, str] = {}

    for family in sorted(families):
        members = families[family]
        genome_counts = Counter(row["genome"] for row in members)
        duplicate_genomes = sorted(g for g, n in genome_counts.items() if n > 1)
        n_hight = sum(row["ecotype"] == "HighT" for row in members)
        n_lowt = sum(row["ecotype"] == "LowT" for row in members)
        n_other = sum(row["ecotype"] == "OtherNP" for row in members)
        missing_cds = 0
        missing_protein = 0
        mismatches = 0
        internal_stop_sequences = 0
        frame_failures = 0
        source_mismatches = 0
        cds_records = []
        protein_records = []
        clean_cds_records = []
        clean_protein_records = []
        clean_ecotypes: list[str] = []

        for row in sorted(members, key=lambda x: x["genome"]):
            locus = row["locus_tag"]
            genome = row["genome"]
            tree_id = tree_id_from_locus(locus)
            previous_tree_id = genome_to_tree_id.setdefault(genome, tree_id)
            if previous_tree_id != tree_id:
                raise ValueError(
                    f"Genome {genome} maps to multiple tree IDs: {previous_tree_id}, {tree_id}"
                )
            cds = ffn.get(locus)
            protein = faa.get(locus)
            if cds is None:
                missing_cds += 1
            if protein is None:
                missing_protein += 1
            status = "not_tested"
            internal_stops = ""
            translated_length = ""
            aa_mismatch_count = ""
            first_aa_mismatch_position = ""
            first_aa_mismatch_pair = ""
            translated_minus_protein_length = ""
            if cds is not None and protein is not None:
                (
                    translated, status, internal_stops, aa_mismatch_count,
                    first_aa_mismatch_position, first_aa_mismatch_pair,
                    translated_minus_protein_length,
                ) = translation_status(cds, protein)
                translated_length = len(translated.rstrip("*"))
                mismatches += status == "mismatch"
                internal_stop_sequences += internal_stops > 0
                frame_failures += len(cds) % 3 != 0
                source_mismatches += (
                    ffn_source.get(locus) != genome or faa_source.get(locus) != genome
                )
                description = f"genome={genome} locus={locus} ecotype={row['ecotype']}"
                cds_records.append((tree_id, description, cds))
                protein_records.append((tree_id, description, protein.rstrip("*")))
                sequence_pass = (
                    status != "mismatch"
                    and internal_stops == 0
                    and len(cds) % 3 == 0
                    and ffn_source.get(locus) == genome
                    and faa_source.get(locus) == genome
                )
                if sequence_pass:
                    clean_cds_records.append((tree_id, description, cds))
                    clean_protein_records.append((tree_id, description, protein.rstrip("*")))
                    clean_ecotypes.append(row["ecotype"])
            else:
                sequence_pass = False
            sequence_rows.append({
                "split_family": family,
                "genome": genome,
                "tree_id": tree_id,
                "ecotype": row["ecotype"],
                "locus_tag": locus,
                "cds_found": int(cds is not None),
                "protein_found": int(protein is not None),
                "cds_length": len(cds) if cds is not None else "",
                "protein_length": len(protein.rstrip("*")) if protein is not None else "",
                "translated_length": translated_length,
                "cds_mod3": len(cds) % 3 if cds is not None else "",
                "internal_stops": internal_stops,
                "translation_status": status,
                "sequence_pass_qc": int(sequence_pass),
                "first_codon": cds[:3] if cds is not None else "",
                "translated_first_aa": translated[0] if cds is not None and protein is not None and translated else "",
                "protein_first_aa": protein[0] if protein else "",
                "aa_mismatch_count": aa_mismatch_count,
                "first_aa_mismatch_position": first_aa_mismatch_position,
                "first_aa_mismatch_pair": first_aa_mismatch_pair,
                "translated_minus_protein_length": translated_minus_protein_length,
                "ffn_source_matches_genome": int(ffn_source.get(locus) == genome),
                "faa_source_matches_genome": int(faa_source.get(locus) == genome),
            })

        stem = safe_name(family)
        with (cds_dir / f"{stem}.fna").open("w", encoding="utf-8") as handle:
            for tree_id, description, seq in cds_records:
                handle.write(f">{tree_id} {description}\n{wrap(seq)}\n")
        with (protein_dir / f"{stem}.faa").open("w", encoding="utf-8") as handle:
            for tree_id, description, seq in protein_records:
                handle.write(f">{tree_id} {description}\n{wrap(seq)}\n")

        with (clean_cds_dir / f"{stem}.fna").open("w", encoding="utf-8") as handle:
            for tree_id, description, seq in clean_cds_records:
                handle.write(f">{tree_id} {description}\n{wrap(seq)}\n")
        with (clean_protein_dir / f"{stem}.faa").open("w", encoding="utf-8") as handle:
            for tree_id, description, seq in clean_protein_records:
                handle.write(f">{tree_id} {description}\n{wrap(seq)}\n")

        clean_counts = Counter(clean_ecotypes)
        eligible = (
            not duplicate_genomes
            and clean_counts["HighT"] >= args.min_hight
            and clean_counts["LowT"] >= args.min_lowt
            and len(clean_cds_records) == len(clean_protein_records)
        )
        if eligible:
            eligible_families.append(family)
        family_rows.append({
            "split_family": family,
            "file_stem": stem,
            "n_sequences": len(members),
            "n_genomes": len(genome_counts),
            "n_HighT": n_hight,
            "n_LowT": n_lowt,
            "n_OtherNP": n_other,
            "n_sequences_after_qc": len(clean_cds_records),
            "n_HighT_after_qc": clean_counts["HighT"],
            "n_LowT_after_qc": clean_counts["LowT"],
            "n_OtherNP_after_qc": clean_counts["OtherNP"],
            "excluded_sequence_count": len(members) - len(clean_cds_records),
            "duplicate_genome_count": len(duplicate_genomes),
            "missing_cds": missing_cds,
            "missing_protein": missing_protein,
            "translation_mismatch_count": mismatches,
            "internal_stop_sequence_count": internal_stop_sequences,
            "frame_failure_count": frame_failures,
            "source_mismatch_count": source_mismatches,
            "paml_input_eligible": int(eligible),
        })

    seq_fields = [
        "split_family", "genome", "tree_id", "ecotype", "locus_tag", "cds_found", "protein_found",
        "cds_length", "protein_length", "translated_length", "cds_mod3", "internal_stops",
        "translation_status", "first_codon", "translated_first_aa", "protein_first_aa",
        "sequence_pass_qc",
        "aa_mismatch_count", "first_aa_mismatch_position", "first_aa_mismatch_pair",
        "translated_minus_protein_length", "ffn_source_matches_genome", "faa_source_matches_genome",
    ]
    fam_fields = [
        "split_family", "file_stem", "n_sequences", "n_genomes", "n_HighT", "n_LowT", "n_OtherNP",
        "n_sequences_after_qc", "n_HighT_after_qc", "n_LowT_after_qc", "n_OtherNP_after_qc",
        "excluded_sequence_count",
        "duplicate_genome_count", "missing_cds", "missing_protein", "translation_mismatch_count",
        "internal_stop_sequence_count", "frame_failure_count", "source_mismatch_count",
        "paml_input_eligible",
    ]
    write_tsv(args.outdir / "sequence_qc.tsv", sequence_rows, seq_fields)
    write_tsv(args.outdir / "family_qc.tsv", family_rows, fam_fields)
    mapping_rows = [
        {"genome": genome, "tree_id": tree_id}
        for genome, tree_id in sorted(genome_to_tree_id.items())
    ]
    write_tsv(args.outdir / "genome_tree_id_map.tsv", mapping_rows, ["genome", "tree_id"])
    (args.outdir / "eligible_families.txt").write_text(
        "\n".join(eligible_families) + ("\n" if eligible_families else ""), encoding="utf-8"
    )
    summary = [
        "# Soft-core sequence extraction summary",
        "",
        f"- Families requested: {len(families)}",
        f"- Member sequences requested: {sum(len(x) for x in families.values())}",
        f"- Exact-run FFN loci indexed: {len(ffn)}",
        f"- Exact-run FAA loci indexed: {len(faa)}",
        f"- Families passing strict extraction/translation QC: {len(eligible_families)}",
        f"- Families requiring review: {len(families) - len(eligible_families)}",
        "",
        "Individual sequences fail QC if CDS/protein retrieval fails, CDS length is not divisible by "
        "three, an internal stop or translation mismatch is present, or the source filename does not "
        "match the member genome. A family remains eligible when it is single-copy and the retained "
        "sequences include HighT >= " + str(args.min_hight) + " and LowT >= " + str(args.min_lowt) + ".",
    ]
    (args.outdir / "SEQUENCE_EXTRACTION_SUMMARY.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary[2:]))


if __name__ == "__main__":
    main()
