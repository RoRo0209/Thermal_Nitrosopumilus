# Phylogeny, OGT and ancestral-state analysis

This module contains the substantive source files recovered for the V4 paper. It is a **partial reconstruction of the workflow**, not a claim that the complete V4 tree/ASR/dating pipeline has been recovered. Files retain their historical names; their figure numbers are not the V4 numbering.

## Workspace and execution

Set `THERMAL_PROJECT_ROOT` to a separate analysis workspace. Place the inputs below in its `1.phylogenetic/` directory. Run scripts from the repository, for example:

```sh
export THERMAL_PROJECT_ROOT=/absolute/path/to/analysis-workspace
Rscript analysis/1.phylogenetic/Final_ASR_clade_analysis.R
Rscript analysis/1.phylogenetic/phylosignal_all108.R
Rscript analysis/1.phylogenetic/Clade_OGT_Explained_Variance.R
Rscript analysis/1.phylogenetic/OGT_Bimodality_Analysis.R
```

In PowerShell use `$env:THERMAL_PROJECT_ROOT = 'D:/analysis-workspace'`. Scripts write outputs into the supplied workspace and may replace outputs with the same names. Input data and generated figures are not included in this repository. No dependency is installed automatically.

## Recovered entry points

| Script | Required inputs | Computation and outputs | V4 relationship / limit |
|---|---|---|---|
| `Final_ASR_clade_analysis.R` | `gtdbtk.ar53.NP.treefile`; `chiplot-tree.xlsx` with `FileName`, `predOGT (C)` / its `openxlsx` name `predOGT.(C)`, and `group` | Reads supplied lowT/highT groups, computes focal-clade OGT comparisons and signal; writes `ASR_clade_assignments_final.csv`, member lists and `Figure1_ASR_clades_final.pdf/png` | OGT and supplied membership processing supporting Fig. 2; **does not infer clade boundaries**. Its focal-clade signal must not be confused with the all-genome analysis. |
| `phylosignal_all108.R` | Same ML tree and `ASR_clade_assignments_final.csv` | Prints Blomberg's K and Pagel's lambda for all-genome ML tree, focal-clade ML tree and chronos tree; K uses 1,000 permutations | Fig. 2c computational source. Select the manuscript-specified tree/subset; the alternatives are not interchangeable. Historical stochastic state was not recorded here. |
| `Clade_OGT_Explained_Variance.R` | `ASR_clade_assignments_final.csv` | Calculates OGT variation explained by supplied clades; prints statistics | Supporting OGT comparison. |
| `OGT_Bimodality_Analysis.R` | `ASR_clade_assignments_final.csv` with `OGT` and `Clade` | Hartigan dip and Gaussian-mixture analysis; `OGT_Bimodality_Distribution.pdf/png`, `OGT_Bimodality_Results.txt` | Supplementary OGT distribution analysis (S5). |
| `OGT_Validation_Analysis_FINAL_CORRECTED.R` | Workspace-root `representative_NP.xlsx`; `cultured_strains_validation_FINAL.csv` | Compares predicted and measured culture OGT, computes correlation, regression, RMSE and MAE, generates validation plots | Measurements/references must be supplied in the input table; static historical strain-result messages were removed from the code. |
| `Generate_FigS1_AB_combined.R` | ML tree, `chiplot-tree.xlsx`, `ASR/NP_nearshore_coastal_abundance.xlsx`, workspace-root `representative_NP.xlsx` | Four-state ARD stochastic mapping, 50 maps, combining median OGT and breadth; node probability output and combined ASR plots | Recovered ASR source, **not established as the exact S4 implementation**. Uses a median OGT threshold and four states; V4 describes fixed 27.5°C thermal states and crown marginal probability >0.80. |
| `Figure1_With_ggtree.R` | `MCMCtree_cleaned.tre` or original `MCMCtree_FigTree.tre`; `chiplot-tree.xlsx`; `TIME2.xlsx` | Draws an already dated tree and time annotations | Related to Fig. 6a/S10. It does not run MCMCtree or recover calibration priors; the exact final composite remains unavailable. |

Tree tip identifiers must match `FileName`; display names and clade membership must come from the supplied metadata. Supply the actual metadata from the target analysis, not labels inferred from R node numbers. These scripts contain no ALE reconciliation.

## Dependencies and gaps

R packages: `ape`, `phytools`, `ggtree`, `ggplot2`, `dplyr`, `openxlsx`, `patchwork`, `scales`, `readxl`, `diptest`, `mclust`, `gridExtra`, `ggrepel` (plus R standard packages). Versions of the original analysis environment were not fully recorded; syntax validation of this release is not numerical reproduction.

Missing or not reconciled: the original marker-selection/alignment/tree-inference commands; TOME 2.0 prediction invocation and trained resources; the exact fixed-threshold V4 binary ASR script and validated membership derivation; Fig. 2a tree/mean-abundance/96-sample heatmap assembly; MCMCtree configuration, calibration and chain diagnostics; exact final Fig. 6a assembly. Do not use the four-state script as an undocumented substitute for binary ASR, and do not report an older K/lambda summary as the result of this release.

The separate four-state SecSSE branch was not selected as a V4 reproduction target. No historical script collection or session dump is included.
