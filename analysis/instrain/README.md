# inStrain profiling and available V4 panel source

Set `THERMAL_PROJECT_ROOT` to a separate workspace containing `instrain/`. This repository contains code only: references, paired reads, sample manifests, profiles, workbooks and result tables are not included.

`run_instrain.sh` maps the listed samples and runs inStrain, writing `bam_instrain/`, `instrain_profile/` and `logs_instrain/`. It requires `fst_26samples.tsv` (sample identifier and clade), `sample_manifest.tsv` (sample identifier followed by R1/R2 paths), and three explicit inputs:

| Variable | Input |
|---|---|
| `THERMAL_AOA_REFERENCE_FASTA` | concatenated matched AOA contigs, originally `all_AOA_contig.fa` |
| `THERMAL_AOA_STB` | `genomes.stb` contig-to-genome mapping |
| `THERMAL_AOA_GENE_FASTA` | matched gene FASTA, originally `all_AOA_contig_trans_gene.fa` |

The source can continue after a missing read pair; audit sample/profile completeness. Shell commands retain legacy unquoted path arguments, so use paths without whitespace. The exact package versions and upstream aggregation/filter provenance require the original environment and metadata.

`out/instrain_boxplot.R` is the available legacy source for V4 Fig. 5c–f. It reads `out/instrain_coverage_dayu1.xlsx`, sheet 1, with `clade`, `nucl_diversity_rarefied`, `SNV_count`, `length`, `pN_pS`, `d_prime_mean`. It computes SNVs per kbp and plots the supplied metrics with the original tests. The upstream profile filtering is not performed by this plotting script.

**The pN/pS panel is provenance-blocked.** The project audit identifies the legacy `pN_pS` field as a nonsynonymous/synonymous count ratio without normalization for site opportunities. The available plotting code does not validate or repair that field. Its inclusion records the source used for the V4 panel, not a claim that the panel is a validated opportunity-corrected pN/pS comparison. A replacement requires original gene-level opportunity counts and an author-approved analysis revision.

`reanalyze_population_metrics.py` is a distinct sensitivity analysis. It excludes the legacy pN/pS column, applies its own coverage/breadth/linkage QC, aggregates sample-by-ecotype metrics, and writes `reanalysis/`. In addition to the workbook fields above it needs sample identifiers, `coverage`, `breadth`, `linked_SNV_count` and the columns checked in its schema guard. It uses summed analysed sequence length for SNV density, length-weighted diversity and linked-pair-weighted D′. These are not the same observational unit or tests as the V4 profile-level figure. Do not silently replace the V4 panel with these results.

Requirements: Bash, Bowtie2, samtools, inStrain; plotting packages `readxl`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `rstatix`; sensitivity Python packages NumPy, pandas >=2.2, SciPy and openpyxl. Source syntax was checked. No profiling, manuscript-value reconciliation or scientific validation is claimed.
