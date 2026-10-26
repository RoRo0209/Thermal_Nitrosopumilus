# Codon models and V4 Figure 6c / Supplementary Figure S12

The numbered scripts retain their historical names. `Figure7c` in filenames means **V4 Figure 6c**. The main panel uses the ratio of two separately fitted foreground omega estimates. It is not itself a direct likelihood-ratio test of equal HighT/LowT omega. Scripts 35–42 perform the separate, direct three-category comparison and convergence follow-up; script 45 provides its supplementary comparison. The later alternative main-panel design was omitted.

## Inputs and execution

Use a separate analysis workspace. Python scripts expose explicit `--help` arguments; Bash wrappers print their required arguments when invoked without them. Pass data paths explicitly. For scripts 16, 28, 29, 34 and 45, pass the analysis workspace as the root argument (`--project-root` for 16). Run scripts 30–32 from that workspace because their stage paths are relative to the working directory. Code can remain in this repository. Linux/WSL and the external command-line tools are required for alignment and codeml execution.

| Stage | Recovered scripts | Required material / next product |
|---|---|---|
| Family selection | 01 | Roary split-family CSV, distinct HighT/LowT genome lists, locus-tag membership crosswalk. The historical `soft85` runs require an explicit `--soft-fraction 0.85`; the source default is 0.90. Confirm the original membership before interpreting identical family counts. |
| Sequences | 05 | Selected membership plus matched CDS/protein FASTA, producing QC FASTAs. |
| Alignments | 06a, optional 06b | MAFFT, PAL2NAL and trimAl; optional IQ-TREE gene trees. The source expects 112 pre-QC families, distinct from the 110 codeml-ready families. |
| Foreground labeling | 06_prepare_family_species_trees | Rooted species tree and verified genome-group lists; pruned per-family trees, within-crown and stem foregrounds. |
| Separate-foreground fitting | 07, 08, 09, 10, 12 | PAML codeml jobs, replicate fits, summaries and convergence audit. Script 12 additionally needs the original candidate integration table via `--integrated-table`; its source integration depends on other unbundled analysis tables. This is an explicit input gap, not an automatically generated step in this release. |
| Summary and annotation | 16, 28, 29, 30, 31, 32 | Historical fitted-summary/audit directories, representative FASTA, KOfam results and BlastKOALA follow-up results. External annotation tools/databases are not bundled. Manual family-level annotation review categories in 32 are retained and are specific to the original 110-family dataset. |
| Direct model comparison | 35, 36, 37, 38 | Reuses foreground trees and alignments. Separate null/alternative three-category fits, multiple starts and audit against the annotated separate-foreground panel. |
| Targeted convergence check | 40, 41, 42 | Follow-up for iscU2; a follow-up script does not establish that its original fits converged. |
| Main/supplementary drawing | 34, 45 | Main separate-foreground panel and supplementary model concordance. Both require the completed and audited tables, including the direct-hit/exclusion tables. |

The historical directories read by 16 are `selection_analysis/downloaded_codeml_results_20260721/selection_analysis/results/09_codeml_summary_soft85/` and `selection_analysis/codeml_analysis_20260721/convergence_review/`. Populate them with the corresponding generated summaries and audit tables, not unrelated runs. Script 29 also expects trimmed proteins in `selection_analysis/workstation_stage_inspect/selection_analysis/results/07_alignments_gene_trees_soft85/protein_trimmed/`; inspect its `aln_dir` before transferring regenerated outputs. These historic path conventions are retained so their origin remains traceable.

Typical downstream calls, after supplying those exact inputs:

```bash
python /path/to/code/analysis/selection_analysis/scripts/16_analyze_plot_dnds_history.py --project-root /path/to/workspace
Rscript /path/to/code/analysis/selection_analysis/scripts/28_plot_panelc_ranked_ratio_only.R /path/to/workspace
Rscript /path/to/code/analysis/selection_analysis/scripts/29_prepare_figure7c_kegg_submission.R /path/to/workspace
# After KOfam and BlastKOALA results, and the direct-model convergence audit:
Rscript /path/to/code/analysis/selection_analysis/scripts/34_plot_panelc_annotated_R.R /path/to/workspace
Rscript /path/to/code/analysis/selection_analysis/scripts/45_plot_supplementary_direct_omega_R.R /path/to/workspace
```

Keep checks for 110 families, successful fits, multiple-start agreement, FDR recalculation and the expected excluded set. A failure is a provenance or analysis discrepancy to resolve, not a reason to remove the check. The release was syntax-checked; the original alignments, model fits, numerical results and rendered panels were not regenerated. Exact tool/database versions and the original candidate integration input remain requirements. Population pN/pS and long-term codon-model dN/dS are different analyses.
