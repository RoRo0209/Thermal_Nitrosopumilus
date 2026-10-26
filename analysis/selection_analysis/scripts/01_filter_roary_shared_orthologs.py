#!/usr/bin/env python3
"""Filter HighT/LowT shared single-copy families from Roary split output.

The script keeps Roary split-family identifiers as the unit for population-
genomic and codon-model analyses, and links those families back to the older
Roary nosplit/ALE family identifiers through shared locus tags.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


ROARY_METADATA_COLUMNS = 14


def read_list(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def split_cell(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in text.split("\t") if item.strip()]


def join_unique(values: list[str]) -> str:
    return ";".join(dict.fromkeys(values))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--roary-csv", type=Path, required=True)
    parser.add_argument("--hight-list", type=Path, required=True)
    parser.add_argument("--lowt-list", type=Path, required=True)
    parser.add_argument("--ale-membership", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--soft-fraction", type=float, default=0.90)
    parser.add_argument("--min-other-for-paml", type=int, default=10)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    roary = pd.read_csv(args.roary_csv, dtype=str, keep_default_na=False)
    genome_columns = list(roary.columns[ROARY_METADATA_COLUMNS:])
    hight = read_list(args.hight_list)
    lowt = read_list(args.lowt_list)

    missing_h = sorted(set(hight) - set(genome_columns))
    missing_l = sorted(set(lowt) - set(genome_columns))
    overlap = sorted(set(hight) & set(lowt))
    if missing_h or missing_l or overlap:
        raise SystemExit(
            "Genome-list validation failed:\n"
            f"  missing HighT: {missing_h}\n"
            f"  missing LowT: {missing_l}\n"
            f"  overlap: {overlap}"
        )

    focal = hight + lowt
    other = [genome for genome in genome_columns if genome not in set(focal)]
    soft_pct = int(round(args.soft_fraction * 100))
    high_cutoff = int(len(hight) * args.soft_fraction + 0.999999)
    low_cutoff = int(len(lowt) * args.soft_fraction + 0.999999)

    family_rows: list[dict[str, object]] = []
    member_rows: list[dict[str, object]] = []
    split_members: dict[str, list[str]] = {}

    for _, row in roary.iterrows():
        family = row["Gene"]
        by_genome = {genome: split_cell(row[genome]) for genome in genome_columns}
        split_members[family] = [gene for genes in by_genome.values() for gene in genes]

        for genome, genes in by_genome.items():
            ecotype = "HighT" if genome in hight else "LowT" if genome in lowt else "OtherNP"
            for gene in genes:
                member_rows.append(
                    {
                        "split_family": family,
                        "genome": genome,
                        "ecotype": ecotype,
                        "locus_tag": gene,
                    }
                )

        def summarize(group: list[str]) -> tuple[int, int, int, int]:
            counts = [len(by_genome[g]) for g in group]
            return (
                sum(count > 0 for count in counts),
                sum(counts),
                max(counts, default=0),
                sum(count > 1 for count in counts),
            )

        hp, hs, hm, hmulti = summarize(hight)
        lp, ls, lm, lmulti = summarize(lowt)
        op, os, om, omulti = summarize(other)
        gp, gs, gm, gmulti = summarize(genome_columns)

        focal_single_copy = hm <= 1 and lm <= 1
        global_single_copy = gm <= 1
        strict = focal_single_copy and hp == len(hight) and lp == len(lowt)
        soft = focal_single_copy and hp >= high_cutoff and lp >= low_cutoff
        soft_global = soft and global_single_copy
        paml_ready = soft_global and op >= args.min_other_for_paml

        family_rows.append(
            {
                "split_family": family,
                "non_unique_gene_name": row.get("Non-unique Gene name", ""),
                "annotation": row.get("Annotation", ""),
                "roary_no_isolates": row.get("No. isolates", ""),
                "roary_no_sequences": row.get("No. sequences", ""),
                "HighT_present": hp,
                "HighT_total_sequences": hs,
                "HighT_max_copy": hm,
                "HighT_multicopy_genomes": hmulti,
                "LowT_present": lp,
                "LowT_total_sequences": ls,
                "LowT_max_copy": lm,
                "LowT_multicopy_genomes": lmulti,
                "OtherNP_present": op,
                "OtherNP_total_sequences": os,
                "OtherNP_max_copy": om,
                "OtherNP_multicopy_genomes": omulti,
                "global_present": gp,
                "global_total_sequences": gs,
                "global_max_copy": gm,
                "global_multicopy_genomes": gmulti,
                "strict_core_focal_singlecopy": strict,
                "soft_core_focal_singlecopy": soft,
                "soft_core_global_singlecopy": soft_global,
                "paml_ready_soft_core": paml_ready,
            }
        )

    families = pd.DataFrame(family_rows)
    members = pd.DataFrame(member_rows)

    ale = pd.read_csv(args.ale_membership, sep="\t", dtype=str, keep_default_na=False)
    required = {"name", "id"}
    if not required.issubset(ale.columns):
        raise SystemExit(f"ALE membership lacks required columns {required}: {list(ale.columns)}")

    locus_to_ale: dict[str, list[str]] = defaultdict(list)
    for family, locus in ale[["name", "id"]].itertuples(index=False):
        if locus:
            locus_to_ale[locus].append(family)

    split_to_parents: dict[str, list[str]] = {}
    parent_to_splits: dict[str, set[str]] = defaultdict(set)
    crosswalk_counts: dict[str, Counter[str]] = {}
    for split_family, loci in split_members.items():
        parent_counts: Counter[str] = Counter()
        for locus in loci:
            parent_counts.update(locus_to_ale.get(locus, []))
        crosswalk_counts[split_family] = parent_counts
        parents = sorted(parent_counts)
        split_to_parents[split_family] = parents
        for parent in parents:
            parent_to_splits[parent].add(split_family)

    crosswalk_rows: list[dict[str, object]] = []
    for split_family, loci in split_members.items():
        parents = split_to_parents[split_family]
        counts = crosswalk_counts[split_family]
        mapped_loci = sum(1 for locus in loci if locus in locus_to_ale)
        if not parents:
            status = "unmapped"
        elif len(parents) > 1:
            status = "mixed_nosplit_parents"
        elif len(parent_to_splits[parents[0]]) == 1:
            status = "one_to_one_family"
        else:
            status = "split_subfamily_of_nosplit_parent"
        crosswalk_rows.append(
            {
                "split_family": split_family,
                "nosplit_ALE_families": ";".join(parents),
                "nosplit_parent_member_counts": ";".join(f"{p}:{counts[p]}" for p in parents),
                "n_split_members": len(loci),
                "n_members_mapped_to_ALE": mapped_loci,
                "member_mapping_fraction": mapped_loci / len(loci) if loci else 0.0,
                "crosswalk_status": status,
            }
        )
    crosswalk = pd.DataFrame(crosswalk_rows)
    families = families.merge(crosswalk, on="split_family", how="left")
    members = members.merge(
        ale[["id", "name", "ko", "kegg"]].rename(
            columns={"id": "locus_tag", "name": "nosplit_ALE_family"}
        ),
        on="locus_tag",
        how="left",
    )

    families = families.sort_values(
        ["soft_core_focal_singlecopy", "HighT_present", "LowT_present", "split_family"],
        ascending=[False, False, False, True],
    )
    families.to_csv(args.outdir / "roary_split_family_qc_all.tsv", sep="\t", index=False)
    members.to_csv(args.outdir / "roary_split_family_members_long.tsv", sep="\t", index=False)
    crosswalk.to_csv(args.outdir / "split_to_nosplit_ALE_crosswalk.tsv", sep="\t", index=False)

    output_filters = {
        "shared_single_copy_strict_focal.tsv": "strict_core_focal_singlecopy",
        f"shared_single_copy_soft{soft_pct}_focal.tsv": "soft_core_focal_singlecopy",
        f"shared_single_copy_soft{soft_pct}_global.tsv": "soft_core_global_singlecopy",
        f"shared_single_copy_soft{soft_pct}_paml_ready.tsv": "paml_ready_soft_core",
    }
    for filename, column in output_filters.items():
        subset = families[families[column]].copy()
        subset.to_csv(args.outdir / filename, sep="\t", index=False)

    soft_families = set(families.loc[families["soft_core_focal_singlecopy"], "split_family"])
    members[members["split_family"].isin(soft_families)].to_csv(
        args.outdir / f"shared_single_copy_soft{soft_pct}_members_long.tsv", sep="\t", index=False
    )

    status_counts = crosswalk["crosswalk_status"].value_counts().to_dict()
    soft = families[families["soft_core_focal_singlecopy"]]
    summary = f"""# Roary split shared-ortholog filtering summary

- Roary families: {len(families):,}
- Roary genomes: {len(genome_columns)}
- HighT genomes matched: {len(hight)}
- LowT genomes matched: {len(lowt)}
- Other Nitrosopumilus genomes: {len(other)}
- Soft-core thresholds: HighT >= {high_cutoff}/{len(hight)}; LowT >= {low_cutoff}/{len(lowt)}
- Strict focal single-copy families: {int(families['strict_core_focal_singlecopy'].sum()):,}
- Soft{soft_pct} focal single-copy families: {int(families['soft_core_focal_singlecopy'].sum()):,}
- Soft{soft_pct} globally single-copy families: {int(families['soft_core_global_singlecopy'].sum()):,}
- PAML-ready soft{soft_pct} families (globally single-copy; other NP >= {args.min_other_for_paml}): {int(families['paml_ready_soft_core'].sum()):,}
- Soft{soft_pct} families with complete locus-to-ALE mapping: {int((soft['member_mapping_fraction'] == 1).sum()):,}/{len(soft):,}

## Split-to-nosplit/ALE crosswalk status

"""
    for status, count in sorted(status_counts.items()):
        summary += f"- {status}: {count:,}\n"
    (args.outdir / "ROARY_SHARED_ORTHOLOG_SUMMARY.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
