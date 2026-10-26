#!/usr/bin/env python3
"""
Analyze SCM1 reciprocal shared-coordinate variant data.

Inputs expected from 01_run_SCM1_mapping_calling.sh:
  results/SCM1_joint_26samples.vcf.gz
  results/SCM1_mapping_QC.tsv
  results/SCM1_amoABC_depth_by_sample.tsv

Main outputs are written under results/ and figures/.
"""

from __future__ import annotations

# Release path adapter: data and outputs live in a separate user workspace.
import os as _release_os
from pathlib import Path as _ReleasePath
if not _release_os.environ.get("THERMAL_PROJECT_ROOT"):
    raise SystemExit("Set THERMAL_PROJECT_ROOT to the analysis workspace")
_release_root = _ReleasePath(_release_os.environ["THERMAL_PROJECT_ROOT"]).expanduser().resolve(strict=True)
_release_os.chdir(_release_root / 'tajima/SCM1_reciprocal_FST')


import csv
import gzip
import math
import os
import random
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str((_release_root / "tajima" / "SCM1_reciprocal_FST") / ".matplotlib_cache"))

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = (_release_root / "tajima" / "SCM1_reciprocal_FST")
RESULTS = BASE / "results"
FIGURES = BASE / "figures"
GROUPS = BASE / "sample_groups_26.tsv"
REF_MANIFEST = BASE / "reference_manifest.tsv"
AMO_LOCI = BASE / "SCM1_amoABC_loci.tsv"
JOINT_VCF = RESULTS / "SCM1_joint_26samples.vcf.gz"
MAPPING_QC = RESULTS / "SCM1_mapping_QC.tsv"
AMO_DEPTH = RESULTS / "SCM1_amoABC_depth_by_sample.tsv"

MIN_DP = int(os.environ.get("MIN_DP", 3))
MIN_CALLABLE_FRACTION = float(os.environ.get("MIN_CALLABLE_FRACTION", 0.80))
WINDOW_SIZE = 10_000
STEP_SIZE = 5_000
N_PERM = 9999
SEED = 20260630

HIGH = "#C85A3A"
LOW = "#2B8C9E"
GREY = "#666666"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 8,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
    }
)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def parse_format(format_keys: list[str], sample_value: str) -> dict[str, str]:
    values = sample_value.split(":")
    return {k: values[i] if i < len(values) else "" for i, k in enumerate(format_keys)}


def parse_number_list(value: str) -> list[float]:
    if not value or value == ".":
        return []
    out = []
    for item in value.split(","):
        try:
            out.append(float(item))
        except ValueError:
            out.append(math.nan)
    return out


def af_from_sample(format_keys: list[str], sample_value: str) -> tuple[float, float]:
    rec = parse_format(format_keys, sample_value)
    dp = math.nan
    if rec.get("DP") not in (None, "", "."):
        try:
            dp = float(rec["DP"])
        except ValueError:
            dp = math.nan

    if "AO" in rec and "RO" in rec:
        ao = parse_number_list(rec.get("AO", ""))
        ro = parse_number_list(rec.get("RO", ""))
        alt = ao[0] if ao else 0.0
        ref = ro[0] if ro else 0.0
        total = ref + alt
        if total > 0:
            return alt / total, total

    if "AD" in rec:
        ad = parse_number_list(rec.get("AD", ""))
        if len(ad) >= 2:
            total = ad[0] + ad[1]
            if total > 0:
                return ad[1] / total, total

    gt = rec.get("GT", "")
    if gt in {"0", "0/0", "0|0"}:
        return 0.0, dp
    if gt in {"1", "0/1", "1/0", "0|1", "1|0", "1/1", "1|1"}:
        if gt in {"1", "1/1", "1|1"}:
            return 1.0, dp
        return 0.5, dp
    return math.nan, dp


def load_vcf_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not JOINT_VCF.exists():
        raise FileNotFoundError(f"Missing joint VCF: {JOINT_VCF}")

    opener = gzip.open if JOINT_VCF.suffix == ".gz" else open
    samples = []
    rows = []
    site_rows = []

    with opener(JOINT_VCF, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, _id, ref, alt = fields[:5]
            if "," in alt:
                continue
            if len(ref) != 1:
                continue
            if alt not in {".", "<*>"} and len(alt) != 1:
                continue
            fmt = fields[8].split(":")
            afs = []
            dps = []
            for sample_value in fields[9:]:
                af, dp = af_from_sample(fmt, sample_value)
                afs.append(af)
                dps.append(dp)
            callable_mask = [not math.isnan(a) and not math.isnan(d) and d >= MIN_DP for a, d in zip(afs, dps)]
            if sum(callable_mask) / len(samples) < MIN_CALLABLE_FRACTION:
                continue
            if np.nanmax(afs) - np.nanmin(afs) <= 0:
                continue
            rows.append(afs)
            site_rows.append(
                {
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "callable_samples": int(sum(callable_mask)),
                    "callable_fraction": float(sum(callable_mask) / len(samples)),
                    "mean_depth": float(np.nanmean(dps)),
                }
            )

    matrix = pd.DataFrame(rows, columns=samples)
    sites = pd.DataFrame(site_rows)
    matrix.insert(0, "site_id", sites["chrom"].astype(str) + ":" + sites["pos"].astype(str))
    return matrix, sites


def sample_distance(matrix: pd.DataFrame) -> pd.DataFrame:
    af = matrix.drop(columns=["site_id"]).T
    samples = list(af.index)
    values = af.to_numpy(dtype=float)
    dist = np.zeros((len(samples), len(samples)), dtype=float)
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            x = values[i]
            y = values[j]
            mask = ~np.isnan(x) & ~np.isnan(y)
            if mask.sum() == 0:
                d = math.nan
            else:
                d = float(np.sqrt(np.mean((x[mask] - y[mask]) ** 2)))
            dist[i, j] = dist[j, i] = d
    return pd.DataFrame(dist, index=samples, columns=samples)


def pcoa(distance: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    d = distance.to_numpy(dtype=float)
    if np.isnan(d).any():
        raise ValueError("Distance matrix contains NaN.")
    n = d.shape[0]
    h = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * h @ (d ** 2) @ h
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = eigvals > 1e-12
    coords = eigvecs[:, positive] * np.sqrt(eigvals[positive])
    explained = eigvals[positive] / eigvals[positive].sum()
    return pd.DataFrame({"sample": distance.index, "PCoA1": coords[:, 0], "PCoA2": coords[:, 1]}), explained


def permanova(distance: pd.DataFrame, groups: pd.Series) -> dict[str, float]:
    labels = groups.to_numpy()
    d = distance.to_numpy(dtype=float)
    n = len(labels)
    d2 = d ** 2
    h = np.eye(n) - np.ones((n, n)) / n
    gmat = -0.5 * h @ d2 @ h
    ss_total = float(np.trace(gmat))

    def f_for(current: np.ndarray) -> tuple[float, float, float]:
        unique = np.unique(current)
        ss_within = 0.0
        for lab in unique:
            idx = np.where(current == lab)[0]
            sub = d2[np.ix_(idx, idx)]
            ss_within += float(sub.sum() / (2.0 * len(idx)))
        ss_between = ss_total - ss_within
        df_between = len(unique) - 1
        df_within = n - len(unique)
        return (ss_between / df_between) / (ss_within / df_within), ss_between, ss_within

    observed_f, ss_between, ss_within = f_for(labels)
    rng = np.random.default_rng(SEED)
    count = 0
    for _ in range(N_PERM):
        f_perm, _, _ = f_for(rng.permutation(labels))
        if f_perm >= observed_f - 1e-12:
            count += 1
    return {
        "method": "one-way PERMANOVA",
        "distance": "SCM1 shared-coordinate allele-frequency Euclidean distance",
        "n_permutations": N_PERM,
        "n_samples": n,
        "n_groups": len(np.unique(labels)),
        "F": observed_f,
        "R2": ss_between / ss_total if ss_total > 0 else math.nan,
        "p_value": (count + 1) / (N_PERM + 1),
    }


def plot_pcoa(coords: pd.DataFrame, explained: np.ndarray, stats: dict[str, float]) -> None:
    fig, ax = plt.subplots(figsize=(3.7, 3.15))
    colors = {"HighT": HIGH, "LowT": LOW}
    for group in ["HighT", "LowT"]:
        sub = coords[coords["group"] == group]
        ax.scatter(
            sub["PCoA1"],
            sub["PCoA2"],
            s=30,
            color=colors[group],
            edgecolor="white",
            linewidth=0.5,
            label=f"{group}-enriched (n={len(sub)})",
        )
    ax.axhline(0, color="#D0D0D0", lw=0.6)
    ax.axvline(0, color="#D0D0D0", lw=0.6)
    ax.set_xlabel(f"PCoA1 ({explained[0] * 100:.1f}%)")
    ax.set_ylabel(f"PCoA2 ({explained[1] * 100:.1f}%)")
    ax.set_title("SCM1-coordinate recruitment PCoA", fontsize=8.5, weight="bold")
    ax.legend(loc="best", fontsize=6.4)
    ax.text(
        0.02,
        0.02,
        f"PERMANOVA: R2={stats['R2']:.3f}, p={stats['p_value']:.4g}",
        transform=ax.transAxes,
        fontsize=6.4,
        color=GREY,
        va="bottom",
    )
    fig.tight_layout()
    fig.savefig(RESULTS / "SCM1_PCoA.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "SCM1_PCoA.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "SCM1_PCoA.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / "SCM1_PCoA.svg", bbox_inches="tight")
    plt.close(fig)


def hudson_like_by_site(matrix: pd.DataFrame, groups: pd.DataFrame, sites: pd.DataFrame) -> pd.DataFrame:
    sample_groups = dict(zip(groups["sample"], groups["group"]))
    sample_cols = [c for c in matrix.columns if c != "site_id"]
    high_cols = [s for s in sample_cols if sample_groups.get(s) == "HighT"]
    low_cols = [s for s in sample_cols if sample_groups.get(s) == "LowT"]
    records = []
    af = matrix.drop(columns=["site_id"])
    for idx, row in af.iterrows():
        h = row[high_cols].astype(float).to_numpy()
        l = row[low_cols].astype(float).to_numpy()
        h = h[~np.isnan(h)]
        l = l[~np.isnan(l)]
        if len(h) < 2 or len(l) < 2:
            continue
        p1 = float(np.mean(h))
        p2 = float(np.mean(l))
        num = (p1 - p2) ** 2 - (p1 * (1 - p1) / (len(h) - 1)) - (p2 * (1 - p2) / (len(l) - 1))
        den = p1 * (1 - p2) + p2 * (1 - p1)
        records.append(
            {
                "site_id": matrix.loc[idx, "site_id"],
                "chrom": sites.loc[idx, "chrom"],
                "pos": int(sites.loc[idx, "pos"]),
                "highT_mean_af": p1,
                "lowT_mean_af": p2,
                "delta_af": p1 - p2,
                "abs_delta_af": abs(p1 - p2),
                "fst_num": num,
                "fst_den": den,
                "site_fst": num / den if den > 0 else math.nan,
                "n_highT_callable": len(h),
                "n_lowT_callable": len(l),
            }
        )
    return pd.DataFrame(records)


def window_fst(site_fst: pd.DataFrame) -> pd.DataFrame:
    records = []
    for chrom, sub in site_fst.groupby("chrom", sort=False):
        max_pos = int(sub["pos"].max())
        for start in range(1, max_pos + 1, STEP_SIZE):
            end = min(start + WINDOW_SIZE - 1, max_pos)
            w = sub[(sub["pos"] >= start) & (sub["pos"] <= end)]
            if w.empty:
                continue
            den = w["fst_den"].sum()
            fst = w["fst_num"].sum() / den if den > 0 else math.nan
            records.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "n_sites": len(w),
                    "Fst": fst,
                    "mean_abs_delta_af": w["abs_delta_af"].mean(),
                }
            )
    out = pd.DataFrame(records)
    if out.empty:
        return out
    threshold = out["Fst"].quantile(0.95)
    out["top5_FST_window"] = out["Fst"] >= threshold
    out["top5_threshold"] = threshold
    return out


def merge_top_windows(windows: pd.DataFrame) -> pd.DataFrame:
    top = windows[windows["top5_FST_window"]].sort_values(["chrom", "start"])
    regions = []
    current = None
    for _, row in top.iterrows():
        if current is None or row["chrom"] != current["chrom"] or row["start"] > current["end"] + STEP_SIZE:
            if current is not None:
                regions.append(current)
            current = {
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "n_windows": 1,
                "max_Fst": float(row["Fst"]),
                "mean_Fst": float(row["Fst"]),
                "mean_abs_delta_af": float(row["mean_abs_delta_af"]),
            }
        else:
            current["end"] = int(max(current["end"], row["end"]))
            current["n_windows"] += 1
            current["max_Fst"] = max(current["max_Fst"], float(row["Fst"]))
            current["mean_Fst"] += float(row["Fst"])
            current["mean_abs_delta_af"] += float(row["mean_abs_delta_af"])
    if current is not None:
        regions.append(current)
    for r in regions:
        r["mean_Fst"] /= r["n_windows"]
        r["mean_abs_delta_af"] /= r["n_windows"]
        r["length_bp"] = r["end"] - r["start"] + 1
    return pd.DataFrame(regions)


def amo_summary(site_fst: pd.DataFrame, windows: pd.DataFrame) -> pd.DataFrame:
    loci = read_tsv(AMO_LOCI)
    depth = read_tsv(AMO_DEPTH) if AMO_DEPTH.exists() else pd.DataFrame()
    rows = []
    threshold = windows["top5_threshold"].dropna().iloc[0] if not windows.empty else math.nan
    for _, locus in loci.iterrows():
        gene = locus["gene"]
        contig = locus["fasta_contig"]
        start = int(locus["start"])
        end = int(locus["end"])
        s = site_fst[(site_fst["chrom"] == contig) & (site_fst["pos"] >= start) & (site_fst["pos"] <= end)]
        w = windows[(windows["chrom"] == contig) & (windows["end"] >= start) & (windows["start"] <= end)]
        d = depth[depth["gene"] == gene] if not depth.empty else pd.DataFrame()
        high_d = d[d["group"] == "HighT"] if not d.empty else pd.DataFrame()
        low_d = d[d["group"] == "LowT"] if not d.empty else pd.DataFrame()
        rows.append(
            {
                "gene": gene,
                "locus_tag": locus["locus_tag"],
                "contig": contig,
                "start": start,
                "end": end,
                "n_callable_snp_sites": len(s),
                "max_abs_delta_af": s["abs_delta_af"].max() if not s.empty else math.nan,
                "mean_abs_delta_af": s["abs_delta_af"].mean() if not s.empty else math.nan,
                "max_site_fst": s["site_fst"].max() if not s.empty else math.nan,
                "overlaps_top5_FST_window": bool((w["Fst"] >= threshold).any()) if not w.empty else False,
                "max_window_Fst": w["Fst"].max() if not w.empty else math.nan,
                "HighT_mean_depth": high_d["mean_depth"].mean() if not high_d.empty else math.nan,
                "HighT_mean_breadth_ge3": high_d["breadth_ge3"].mean() if not high_d.empty and "breadth_ge3" in high_d else math.nan,
                "LowT_mean_depth": low_d["mean_depth"].mean() if not low_d.empty else math.nan,
                "LowT_mean_breadth_ge3": low_d["breadth_ge3"].mean() if not low_d.empty and "breadth_ge3" in low_d else math.nan,
            }
        )
    return pd.DataFrame(rows)


def write_report(
    groups: pd.DataFrame,
    qc: pd.DataFrame,
    matrix: pd.DataFrame,
    coords: pd.DataFrame,
    explained: np.ndarray,
    permanova_stats: dict[str, float],
    windows: pd.DataFrame,
    regions: pd.DataFrame,
    amo: pd.DataFrame,
) -> None:
    spot_summary = Path("../../results/fst_summary.txt")
    spot_line = "SPOT01 genome-wide Hudson's FST was not read in this run."
    if spot_summary.exists():
        text = spot_summary.read_text(encoding="utf-8", errors="replace")
        for line in text.splitlines():
            if "Genome-wide Hudson" in line:
                spot_line = line
                break

    n_high = int((groups["group"] == "HighT").sum())
    n_low = int((groups["group"] == "LowT").sum())
    n_pass = int((qc["passes_original_thresholds"] == "yes").sum()) if "passes_original_thresholds" in qc else 0
    genome_fst = windows["Fst"].mean() if not windows.empty else math.nan
    top_threshold = windows["top5_threshold"].dropna().iloc[0] if not windows.empty else math.nan

    lines = [
        "# SCM1 reciprocal FST validation report",
        "",
        "## Aim",
        "",
        "All 26 metagenomic recruitments used in the SPOT01-coordinate analysis were reciprocally mapped to the HighT representative genome SCM1. This checks whether the SPOT01-based structure is robust to the choice of shared reference coordinate.",
        "",
        "## Reference and samples",
        "",
        "- Shared coordinate: SCM1 HighT representative genome.",
        "- Fasta contig used in the VCF: `NC_010085.1`.",
        "- GFF seqid in the local annotation file: `gnl|X|ANJPDMEM_1`; this differs from the fasta contig name and is recorded in `reference_manifest.tsv`.",
        f"- Samples: {n_high} HighT-enriched and {n_low} LowT-enriched recruitments; no re-filtering of sample identity was performed.",
        "",
        "## Mapping QC",
        "",
        f"- Samples passing the original depth/breadth thresholds after SCM1 mapping: {n_pass}/26.",
        f"- Mean SCM1 depth across samples: {qc['mean_depth'].mean():.3f}.",
        f"- Mean SCM1 breadth >=3x across samples: {qc.filter(like='breadth').iloc[:,0].mean():.3f}.",
        "",
        "## Allele-frequency matrix and PCoA",
        "",
        f"- Callable polymorphic SCM1 sites retained: {matrix.shape[0]}.",
        f"- PCoA1 variance: {explained[0] * 100:.2f}%.",
        f"- PCoA2 variance: {explained[1] * 100:.2f}%.",
        f"- PERMANOVA: R2 = {permanova_stats['R2']:.4f}, p = {permanova_stats['p_value']:.5f}.",
        "",
        "## FST windows",
        "",
        f"- Window size: {WINDOW_SIZE} bp; step size: {STEP_SIZE} bp.",
        f"- Number of windows: {len(windows)}.",
        f"- Mean window FST: {genome_fst:.4f}.",
        f"- Top 5% FST threshold: {top_threshold:.4f}.",
        f"- Candidate high-FST regions after merging adjacent top windows: {len(regions)}.",
        "",
        "## amoABC",
        "",
        "SCM1 amoABC loci were assigned by locus-tag homology because the local Prokka GFF labels these CDS as hypothetical proteins.",
    ]
    if not amo.empty:
        for _, r in amo.iterrows():
            lines.append(
                f"- {r['gene']} ({r['locus_tag']}): callable SNP sites = {r['n_callable_snp_sites']}, "
                f"max |delta AF| = {r['max_abs_delta_af']:.4g}, max site FST = {r['max_site_fst']:.4g}, "
                f"overlaps top5 FST window = {r['overlaps_top5_FST_window']}."
            )
    lines.extend(
        [
            "",
            "## Comparison with SPOT01",
            "",
            f"- SPOT01 reference result: {spot_line}",
            "- Interpret SCM1 and SPOT01 jointly. Concordant PCoA/FST/amoABC signals would support a reference-robust ecotype differentiation signal; weaker SCM1 coverage or reduced LowT mapping would indicate potential HighT-reference bias and should be reported as a limitation.",
            "",
            "## Suggested use",
            "",
            "Use the SCM1 reciprocal result as supplementary validation. In the main text, keep the claim conservative: reciprocal mapping to SCM1 was used to test whether SPOT01-coordinate allele-frequency differentiation was driven by reference choice.",
        ]
    )
    (RESULTS / "SCM1_reciprocal_validation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    groups = read_tsv(GROUPS)
    qc = read_tsv(MAPPING_QC)

    matrix, sites = load_vcf_matrix()
    matrix.to_csv(RESULTS / "SCM1_allele_frequency_matrix.tsv", sep="\t", index=False)
    sites.to_csv(RESULTS / "SCM1_callable_sites.tsv", sep="\t", index=False)

    dist = sample_distance(matrix)
    dist.to_csv(RESULTS / "SCM1_pairwise_AF_distance.tsv", sep="\t")
    coords, explained = pcoa(dist)
    coords = coords.merge(groups, on="sample", how="left")
    coords.to_csv(RESULTS / "SCM1_PCoA_coordinates.tsv", sep="\t", index=False)
    stats = permanova(dist, coords["group"])
    pd.DataFrame([stats]).to_csv(RESULTS / "SCM1_PERMANOVA_report.txt", sep="\t", index=False)
    plot_pcoa(coords, explained, stats)

    site_fst = hudson_like_by_site(matrix, groups, sites)
    site_fst.to_csv(RESULTS / "SCM1_site_FST_deltaAF.tsv", sep="\t", index=False)
    windows = window_fst(site_fst)
    windows.to_csv(RESULTS / "SCM1_FST_windows.tsv", sep="\t", index=False)
    regions = merge_top_windows(windows)
    regions.to_csv(RESULTS / "SCM1_top5_FST_regions.tsv", sep="\t", index=False)
    amo = amo_summary(site_fst, windows)
    amo.to_csv(RESULTS / "SCM1_amoABC_summary.tsv", sep="\t", index=False)

    write_report(groups, qc, matrix, coords, explained, stats, windows, regions, amo)
    print("SCM1 reciprocal analysis complete.")
    print(f"Callable polymorphic sites: {matrix.shape[0]}")
    print(f"PERMANOVA R2={stats['R2']:.4f}, p={stats['p_value']:.5f}")


if __name__ == "__main__":
    main()
