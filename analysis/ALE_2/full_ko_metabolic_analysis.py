from pathlib import Path
import math
import re

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr


import os
import json
if not os.environ.get("THERMAL_PROJECT_ROOT") or not os.environ.get("THERMAL_GENOME_GROUPS_JSON"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT and THERMAL_GENOME_GROUPS_JSON (verified 108-genome cohort)")
ROOT = Path(os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True) / "ALE_2"
OUT = ROOT / "full_ko_analysis"
OUT.mkdir(exist_ok=True)

MATRIX = ROOT / "genome_x_family_copy_number_kegg_annotated.tsv"
MEMBERSHIP = ROOT / "gene_family_membership.tsv"
META = ROOT / "tree-gen-name-gff.xlsx"

# The representative-genome cohort is distinct from the full ALE crown.
groups = json.loads(Path(os.environ["THERMAL_GENOME_GROUPS_JSON"]).read_text(encoding="utf-8"))
LOWT, HIGHT = groups["LowT"], groups["HighT"]
EXCLUDED_FROM_108 = set(groups["excluded_from_108"])
if len(LOWT) != 23 or len(HIGHT) != 11 or len(set(LOWT + HIGHT)) != 34:
    raise ValueError("Expected unique, disjoint V4 representative-genome groups: LowT=23, HighT=11")
if EXCLUDED_FROM_108.intersection(LOWT + HIGHT):
    raise ValueError("Excluded genomes cannot belong to focal groups")


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1)
    return out


def suspicious(label):
    x = str(label).lower()
    terms = ["zbtb", "ankfy", "ing5", "ino80", "memo1", "mmp20", "eukary"]
    return any(t in x for t in terms)


mat = pd.read_csv(MATRIX, sep="\t").set_index("genome")
mem = pd.read_csv(MEMBERSHIP, sep="\t", dtype=str).fillna("")
mem = mem[mem["name"].isin(mat.columns)].copy()
mem["ko"] = mem["ko"].str.strip()
mem = mem[mem["ko"].str.fullmatch(r"K\d{5}")].copy()

# One family can contain many genes, but should map to one KO for this analysis.
family_ko = mem.groupby("name").agg(
    ko=("ko", lambda x: sorted(set(x))[0] if len(set(x)) else ""),
    ko_n=("ko", lambda x: len(set(x))),
    annotation=("kegg", lambda x: " | ".join(sorted(set(v for v in x if v)))[:1000]),
).reset_index()
family_ko = family_ko[family_ko.ko_n == 1].copy()

ko_to_families = family_ko.groupby("ko")["name"].apply(list)
ko_matrix = pd.DataFrame(
    {ko: mat[families].sum(axis=1) for ko, families in ko_to_families.items()},
    index=mat.index,
)

annotation = family_ko.groupby("ko")["annotation"].apply(
    lambda x: " | ".join(sorted(set(v for v in x if v)))[:2000]
).to_dict()

missing_low = sorted(set(LOWT) - set(ko_matrix.index))
missing_high = sorted(set(HIGHT) - set(ko_matrix.index))
if missing_low or missing_high:
    raise RuntimeError(f"Missing focal genomes: LowT={missing_low}, HighT={missing_high}")

meta = pd.read_excel(META, sheet_name="Sheet1")
np_genomes = set(meta.loc[meta["type"].eq("g__Nitrosopumilus"), "gff_name"].dropna().astype(str))
np_genomes &= set(ko_matrix.index)
np_genomes -= EXCLUDED_FROM_108
if len(np_genomes) != 108 or not set(LOWT + HIGHT) <= np_genomes:
    raise ValueError("Metadata, matrix and cohort manifest must define the same 108 genomes")
other_np = sorted(np_genomes - set(LOWT) - set(HIGHT))

rows = []
for ko in ko_matrix.columns:
    h = ko_matrix.loc[HIGHT, ko]
    l = ko_matrix.loc[LOWT, ko]
    hp = int((h > 0).sum())
    lp = int((l > 0).sum())
    table = [[hp, len(HIGHT) - hp], [lp, len(LOWT) - lp]]
    odds, p = fisher_exact(table, alternative="two-sided")
    try:
        _, p_cn = mannwhitneyu(h, l, alternative="two-sided")
    except ValueError:
        p_cn = 1.0
    oprev = float((ko_matrix.loc[other_np, ko] > 0).mean()) if other_np else math.nan
    rows.append({
        "KO": ko,
        "annotation": annotation.get(ko, ""),
        "n_families": len(ko_to_families[ko]),
        "HighT_present": hp,
        "LowT_present": lp,
        "HighT_prevalence": hp / len(HIGHT),
        "LowT_prevalence": lp / len(LOWT),
        "OtherNP_prevalence": oprev,
        "delta_prevalence_H_minus_L": hp / len(HIGHT) - lp / len(LOWT),
        "HighT_mean_copy": float(h.mean()),
        "LowT_mean_copy": float(l.mean()),
        "OtherNP_mean_copy": float(ko_matrix.loc[other_np, ko].mean()) if other_np else math.nan,
        "delta_mean_copy_H_minus_L": float(h.mean() - l.mean()),
        "fisher_odds_ratio": odds,
        "fisher_p": p,
        "copy_number_mw_p": p_cn,
    })

res = pd.DataFrame(rows)
res["fisher_q_bh_all_KOs"] = bh(res.fisher_p)
res["copy_number_q_bh_all_KOs"] = bh(res.copy_number_mw_p)
res["direction"] = np.where(res.delta_prevalence_H_minus_L > 0, "HighT-biased",
                     np.where(res.delta_prevalence_H_minus_L < 0, "LowT-biased", "equal"))
res["effect_abs_delta_ge_0.35"] = res.delta_prevalence_H_minus_L.abs() >= 0.35
res["candidate_q05_delta035"] = (
    (res["fisher_q_bh_all_KOs"] < 0.05)
    & res["effect_abs_delta_ge_0.35"]
)
res["suspicious_annotation"] = res.annotation.map(suspicious)

def context_class(r):
    lo, hi, other = r.LowT_prevalence, r.HighT_prevalence, r.OtherNP_prevalence
    if not np.isfinite(other) or hi == lo:
        return "unresolved"
    pos = (other - lo) / (hi - lo)
    if r.direction == "HighT-biased":
        if pos >= 2/3: return "primarily_LowT_reduction"
        if pos <= 1/3: return "primarily_HighT_expansion"
    elif r.direction == "LowT-biased":
        if pos >= 2/3: return "primarily_LowT_expansion"
        if pos <= 1/3: return "primarily_HighT_reduction"
    return "bidirectional_or_intermediate"

res["genus_context"] = res.apply(context_class, axis=1)
res = res.sort_values(["candidate_q05_delta035", "fisher_q_bh_all_KOs", "delta_prevalence_H_minus_L"], ascending=[False, True, False])
res.to_csv(OUT / "all_KO_HighT_LowT_statistics.tsv", sep="\t", index=False)
ko_matrix.to_csv(OUT / "genome_by_KO_copy_number.tsv", sep="\t")

cand = res[res.candidate_q05_delta035].copy()
cand.to_csv(OUT / "significant_KO_candidates_q05_delta035.tsv", sep="\t", index=False)

# Family-level analysis used by the original Fig. 5. This is intentionally
# separate from the KO-level analysis because several paralogous families can
# share one KO label.
family_annotation = mem.groupby("name").agg(
    KO=("ko", lambda x: ",".join(sorted(set(x)))),
    annotation=("kegg", lambda x: " | ".join(sorted(set(v for v in x if v)))[:2000]),
).reindex(mat.columns).fillna("")
family_rows = []
for family in mat.columns:
    h = mat.loc[HIGHT, family]
    l = mat.loc[LOWT, family]
    hp, lp = int((h > 0).sum()), int((l > 0).sum())
    odds, p = fisher_exact(
        [[hp, len(HIGHT) - hp], [lp, len(LOWT) - lp]],
        alternative="two-sided",
    )
    family_rows.append({
        "family": family,
        "KO": family_annotation.loc[family, "KO"],
        "annotation": family_annotation.loc[family, "annotation"],
        "HighT_present": hp,
        "LowT_present": lp,
        "HighT_prevalence": hp / len(HIGHT),
        "LowT_prevalence": lp / len(LOWT),
        "OtherNP_prevalence": float((mat.loc[other_np, family] > 0).mean()),
        "delta_prevalence_H_minus_L": hp / len(HIGHT) - lp / len(LOWT),
        "HighT_mean_copy": h.mean(),
        "LowT_mean_copy": l.mean(),
        "OtherNP_mean_copy": float(mat.loc[other_np, family].mean()),
        "fisher_odds_ratio": odds,
        "fisher_p": p,
    })
family_res = pd.DataFrame(family_rows)
family_res["fisher_q_bh_all_families"] = bh(family_res.fisher_p)
family_res["direction"] = np.where(
    family_res.delta_prevalence_H_minus_L > 0,
    "HighT-biased",
    np.where(family_res.delta_prevalence_H_minus_L < 0, "LowT-biased", "equal"),
)
family_res["candidate_q05_delta035"] = (
    (family_res.fisher_q_bh_all_families < 0.05)
    & (family_res.delta_prevalence_H_minus_L.abs() >= 0.35)
)
family_res["genus_context"] = family_res.apply(context_class, axis=1)
family_res = family_res.sort_values(
    ["candidate_q05_delta035", "fisher_q_bh_all_families"],
    ascending=[False, True],
)
family_res.to_csv(OUT / "all_gene_family_HighT_LowT_statistics.tsv", sep="\t", index=False)
family_res[family_res.candidate_q05_delta035].to_csv(
    OUT / "significant_gene_family_candidates_q05_delta035.tsv", sep="\t", index=False
)

# Genome-level annotation breadth and size-confounding diagnostics.
breadth = pd.DataFrame({
    "genome": ko_matrix.index,
    "annotated_KO_richness": (ko_matrix > 0).sum(axis=1),
    "annotated_KO_total_copies": ko_matrix.sum(axis=1),
}).set_index("genome")
breadth["group"] = "Non-Nitrosopumilus"
breadth.loc[other_np, "group"] = "Other Nitrosopumilus"
breadth.loc[LOWT, "group"] = "LowT"
breadth.loc[HIGHT, "group"] = "HighT"
breadth.to_csv(OUT / "genome_KO_breadth.tsv", sep="\t")

summary = []
summary.append(f"Input genomes: {len(ko_matrix)}; Nitrosopumilus in metadata: {len(np_genomes)}")
summary.append(f"Groups: HighT={len(HIGHT)}, LowT={len(LOWT)}, other Nitrosopumilus={len(other_np)}")
summary.append(f"Unique KOs tested after family-to-KO aggregation: {len(res)}")
summary.append(f"FDR q<0.05: HighT-biased={len(res[(res.fisher_q_bh_all_KOs<.05)&(res.direction=='HighT-biased')])}, LowT-biased={len(res[(res.fisher_q_bh_all_KOs<.05)&(res.direction=='LowT-biased')])}")
summary.append(f"FDR q<0.05 and |delta prevalence|>=0.35: HighT-biased={len(cand[cand.direction=='HighT-biased'])}, LowT-biased={len(cand[cand.direction=='LowT-biased'])}")
family_cand = family_res[family_res.candidate_q05_delta035]
summary.append(
    "Gene-family FDR q<0.05 and |delta prevalence|>=0.35: "
    f"HighT-biased={len(family_cand[family_cand.direction=='HighT-biased'])}, "
    f"LowT-biased={len(family_cand[family_cand.direction=='LowT-biased'])}"
)
summary.append("Genus-context counts among effect-filtered FDR candidates: " + str(cand.genus_context.value_counts().to_dict()))
summary.append("Suspicious annotation candidates: " + str(int(cand.suspicious_annotation.sum())))
for g in ["HighT", "LowT", "Other Nitrosopumilus"]:
    x = breadth[breadth.group == g]
    summary.append(f"{g} annotated KO richness median={x.annotated_KO_richness.median():.1f}, range={x.annotated_KO_richness.min()}-{x.annotated_KO_richness.max()}, n={len(x)}")
(OUT / "analysis_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))
