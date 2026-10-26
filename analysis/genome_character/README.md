# Genome traits, functional annotation summaries and amino-acid PCA

These are recovered analysis sources associated with V4 Fig. 4a–d. Raw annotations, genome metadata, sequences and result tables are supplied separately. Set `THERMAL_PROJECT_ROOT` to an existing analysis workspace mirroring the paths below; scripts write there, not into this code repository.

| Entry point | Inputs relative to the workspace | Outputs / role |
|---|---|---|
| `two_group/boxplot-7.R` | `genome_character/two_group/two_group.xlsx` | Trait comparison panels (`two_group_equalwidth_fixedorder_with_mu.pdf/png`) covering OGT, GC, genome size, coding density, functional counts and gRodon growth rate; related to Fig. 4a and aggregate functional traits. |
| `two_group/combined_barplots.R` | Same `two_group.xlsx`; `two_group/overview.txt`, `tcdb_blast`, `merops_blast`, `my_mer_families.txt` | CAZyme, TCDB class and MEROPS family counts; Wilcoxon tests with BH adjustment; `combined_family_barplots.pdf/png` (Fig. 4b analytical source). |
| `PCA-AA/AA_PCA_by_group_final_contrib.R` | Explicit environment variables `THERMAL_AA_USAGE` and `THERMAL_AA_METADATA` pointing to the AA-usage TSV and group-metadata XLSX | Centered, scaled amino-acid PCA, contribution-colored loading arrows and `AA_PCA_biplot_contrib.pdf/png` (Fig. 4d analytical source). |
| `summarize_functional_counts_all_np.py` | Workspace-root `representative_NP_Supplementary.xlsx`; `genome_character/overview.txt`, `tcdb_blast`, `merops_blast` | All-NP functional count matrix, clade summaries and pairwise statistics; outputs start `functional_gene_counts_all_np_`. |
| `summarize_cog.py` | Command-line directory containing `.emapper.annotations` files | `cog_summary_for_YAS.tsv` and `yas_s_strategy_indicators.tsv` in `genome_character/`; annotation counting input to trait analysis, **not the final Y–A–S radar calculation**. |

## Input schema

`two_group.xlsx` requires `FileName`, `group` (`lowT` or `highT`), `predOGT`, `GC-content`, `Genome_Size` in bp, `Coding_Density`, `dbcan`, `merops`, `tcdb`, and growth information (`mu_per_hour` or a supported doubling-time column). Coding density at 0–1 scale is converted to percent by the original code. The final paper's grouping must be supplied, not inferred here.

dbCAN `overview.txt` is a data table with `Gene ID`, `#ofTools` and annotation columns; genes require support from at least two tools. Gene IDs use a genome prefix followed by `__`. BLAST tables are tab-delimited with standard 12 columns plus optional `qcovs`; filters preserve the original e-value ≤1e−5 and identity ≥30%. `my_mer_families.txt` supplies MEROPS family mapping.

AA usage requires `Genome ID` plus numerical amino-acid columns; metadata requires `FileName` and `group`. The source's labels are `cold-wide`, `warm-wide`, `WCA`, and `Other`; the correspondence to the manuscript groups must be established in metadata. Do not silently rename or infer group membership. Centering/scaling and selection of 12 loading arrows are preserved.

The all-NP summary reads its detailed master schema directly in code (including genome accession fields and `Clade`). Missing annotations are distinguished from annotated genomes with zero qualifying hits.

## Running

```sh
export THERMAL_PROJECT_ROOT=/absolute/path/to/analysis-workspace
Rscript analysis/genome_character/two_group/boxplot-7.R
Rscript analysis/genome_character/two_group/combined_barplots.R
export THERMAL_AA_USAGE=/absolute/path/to/comparem-aa_usage_group_final.txt
export THERMAL_AA_METADATA=/absolute/path/to/NP_WCA_genome.xlsx
Rscript analysis/genome_character/PCA-AA/AA_PCA_by_group_final_contrib.R
python analysis/genome_character/summarize_functional_counts_all_np.py
python analysis/genome_character/summarize_cog.py /absolute/path/to/eggnog_results
```

R dependencies: `readxl`, `readr`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, `patchwork`. Python: `pandas`, `openpyxl`, plus the standard library. Install these before execution; scripts do not install packages. Historical package/database versions and upstream annotation command lines are not fully recoverable from these sources. The trait plot now stops when neither a finite growth rate nor a positive doubling time is supplied, instead of drawing an empty growth panel.

## Remaining gaps

The exact eight-trait, overall-genome-mean-scaled Y–A–S radar script for Fig. 4c was not found. The COG count script is only an upstream component; legacy suggestions to use z-scores were removed because they do not implement the V4 method. Full upstream gRodon, CompareM, dbCAN, MEROPS/TCDB and eggNOG runs need their original database versions and parameters. The final panel assembly may differ from these standalone plots. Accessory-gene/HSP20/nadM analyses (Fig. 4e/f) are handled separately under `analysis/ALE_2/` if available.
