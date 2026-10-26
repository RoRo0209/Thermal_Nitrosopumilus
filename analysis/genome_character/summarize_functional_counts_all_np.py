import os
from pathlib import Path
_thermal_root = os.environ.get("THERMAL_PROJECT_ROOT", "")
if not _thermal_root or not Path(_thermal_root).is_dir():
    raise SystemExit("Set THERMAL_PROJECT_ROOT to an existing, separate analysis workspace.")
THERMAL_PROJECT_ROOT = Path(_thermal_root).resolve()


import pandas as pd


BASE = THERMAL_PROJECT_ROOT / 'genome_character'
ROOT = BASE.parent
MASTER = ROOT / "representative_NP_Supplementary.xlsx"

OUT_COUNTS = BASE / "functional_gene_counts_all_np_processed.csv"
OUT_SUMMARY = BASE / "functional_gene_counts_all_np_clade_summary.tsv"
OUT_PAIRWISE = BASE / "functional_gene_counts_all_np_pairwise.tsv"
OUT_TEXT = BASE / "functional_gene_counts_all_np_summary.md"


def genome_from_gene(series: pd.Series) -> pd.Series:
    return series.astype(str).str.split("__", n=1).str[0]


def summarize_cazyme() -> tuple[pd.Series, set[str]]:
    overview = pd.read_csv(BASE / "overview.txt", sep="\t")
    overview["Genome_ID"] = genome_from_gene(overview["Gene ID"])
    source_genomes = set(overview["Genome_ID"])
    supported = overview[pd.to_numeric(overview["#ofTools"], errors="coerce") >= 2]
    counts = supported.groupby("Genome_ID")["Gene ID"].nunique()
    return counts, source_genomes


def summarize_blast(filename: str) -> tuple[pd.Series, set[str]]:
    columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qcovs",
    ]
    blast = pd.read_csv(BASE / filename, sep="\t", header=None, names=columns)
    blast["Genome_ID"] = genome_from_gene(blast["qseqid"])
    source_genomes = set(blast["Genome_ID"])
    filtered = blast[
        (pd.to_numeric(blast["evalue"], errors="coerce") <= 1e-5)
        & (pd.to_numeric(blast["pident"], errors="coerce") >= 30)
    ].copy()
    filtered = filtered.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
    counts = filtered.groupby("Genome_ID")["qseqid"].nunique()
    return counts, source_genomes


def add_count(table: pd.DataFrame, column: str, counts: pd.Series, sources: set[str]) -> None:
    present = table["Genome_ID"].isin(sources)
    table[f"{column}_source_status"] = present.map({True: "annotated", False: "not annotated"})
    table[column] = table["Genome_ID"].map(counts)
    table.loc[present & table[column].isna(), column] = 0
    table[column] = table[column].astype("Int64")


def cliffs_delta(x, y) -> float:
    vals = [(a > b) - (a < b) for a in x for b in y]
    return sum(vals) / len(vals)


def mann_whitney_u_p_approx(x, y) -> float:
    # Normal approximation with tie correction, two-sided.
    import math

    x = list(map(float, x))
    y = list(map(float, y))
    n1, n2 = len(x), len(y)
    vals = pd.Series(x + y)
    ranks = vals.rank(method="average").to_numpy()
    r1 = ranks[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2
    u2 = n1 * n2 - u1
    u = min(u1, u2)
    n = n1 + n2
    counts = vals.value_counts().to_numpy()
    tie_term = sum(c**3 - c for c in counts)
    var = n1 * n2 / 12 * ((n + 1) - tie_term / (n * (n - 1)))
    mean = n1 * n2 / 2
    z = (u - mean + 0.5) / math.sqrt(var)
    return math.erfc(abs(z) / math.sqrt(2))


def main() -> None:
    master = pd.read_excel(MASTER)
    table = master[
        [
            "Genome_ID",
            "GTDB_Accession_ID",
            "NCBI_Accession_ID",
            "Strain_Name",
            "Clade",
            "Predicted_OGT_C",
            "Genome_Size",
            "GC_Content",
            "Genome_Completeness",
            "Genome_Contamination",
        ]
    ].copy()

    cazyme_counts, cazyme_sources = summarize_cazyme()
    peptidase_counts, peptidase_sources = summarize_blast("merops_blast")
    transporter_counts, transporter_sources = summarize_blast("tcdb_blast")

    add_count(table, "CAZymes_count", cazyme_counts, cazyme_sources)
    add_count(table, "Peptidases_count", peptidase_counts, peptidase_sources)
    add_count(table, "Transporters_count", transporter_counts, transporter_sources)

    source_cols = [
        "CAZymes_count_source_status",
        "Peptidases_count_source_status",
        "Transporters_count_source_status",
    ]
    table["n_sources_annotated"] = (table[source_cols] == "annotated").sum(axis=1)
    table["Annotation_coverage"] = table["n_sources_annotated"].map(
        {3: "complete: 3/3 sources", 2: "partial: 2/3 sources", 1: "partial: 1/3 sources", 0: "not annotated"}
    )
    table.to_csv(OUT_COUNTS, index=False)

    features = ["CAZymes_count", "Peptidases_count", "Transporters_count"]
    summary_rows = []
    for feature in features:
        for clade, sub in table.groupby("Clade"):
            values = pd.to_numeric(sub[feature], errors="coerce").dropna()
            summary_rows.append(
                {
                    "feature": feature,
                    "clade": clade,
                    "n_total": len(sub),
                    "n_annotated": len(values),
                    "median": values.median() if len(values) else pd.NA,
                    "mean": values.mean() if len(values) else pd.NA,
                    "min": values.min() if len(values) else pd.NA,
                    "max": values.max() if len(values) else pd.NA,
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    pair_rows = []
    pairs = [("highT", "lowT"), ("highT", "other_NP"), ("lowT", "other_NP")]
    for feature in features:
        for a, b in pairs:
            av = pd.to_numeric(table.loc[table["Clade"] == a, feature], errors="coerce").dropna()
            bv = pd.to_numeric(table.loc[table["Clade"] == b, feature], errors="coerce").dropna()
            if len(av) and len(bv):
                pair_rows.append(
                    {
                        "feature": feature,
                        "comparison": f"{a}_vs_{b}",
                        f"{a}_n": len(av),
                        f"{b}_n": len(bv),
                        f"{a}_median": av.median(),
                        f"{b}_median": bv.median(),
                        f"{a}_mean": av.mean(),
                        f"{b}_mean": bv.mean(),
                        "mann_whitney_p_approx": mann_whitney_u_p_approx(av, bv),
                        f"cliffs_delta_{a}_minus_{b}": cliffs_delta(av.tolist(), bv.tolist()),
                    }
                )
    pairwise = pd.DataFrame(pair_rows)
    pairwise.to_csv(OUT_PAIRWISE, sep="\t", index=False)

    cov = table["Annotation_coverage"].value_counts()
    lines = [
        "# All-NP functional count summary",
        "",
        "## Annotation coverage",
        "",
        f"- Master NP genomes: {len(table)}.",
        f"- Complete local MEROPS/dbCAN/TCDB annotation: {int(cov.get('complete: 3/3 sources', 0))}.",
        f"- Complete annotation by clade: {table[table['Annotation_coverage'] == 'complete: 3/3 sources']['Clade'].value_counts().to_dict()}.",
        "",
        "## Clade medians",
        "",
    ]
    for feature in features:
        sub = summary[summary["feature"] == feature]
        med = {r["clade"]: r["median"] for _, r in sub.iterrows()}
        lines.append(
            f"- {feature}: highT {med.get('highT')}, lowT {med.get('lowT')}, other_NP {med.get('other_NP')}."
        )
    lines.extend(["", "## Pairwise tests", ""])
    for _, r in pairwise.iterrows():
        lines.append(
            f"- {r['feature']} {r['comparison']}: p ≈ {r['mann_whitney_p_approx']:.3g}."
        )
    OUT_TEXT.write_text("\n".join(lines), encoding="utf-8")

    print(OUT_TEXT)


if __name__ == "__main__":
    main()
