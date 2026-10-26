# Source packaging changes — 2026-10-26

The target is the V4 manuscript and supplementary figures. Original research files were preserved. The release contains selected available sources rather than a copy of the entire working project.

- Replaced machine-specific paths and private R-library injections with explicit arguments or `THERMAL_PROJECT_ROOT` and documented environmental-data variables. Inputs/outputs belong in a separate workspace.
- Preserved substantive source calculations, model definitions, filtering choices and study-specific consistency checks. Legacy script names were retained; the figure map supplies V4 numbering.
- Added same-run ALE clade matching and per-family event extraction, using actual numbered tree/terminal records and explicit verified crown membership. Added negative tests for node/run mismatches, missing/duplicate branches and invalid group clades. No assumed numeric crown labels bypass these checks.
- Externalized selected reference-region labels and genomic group lists into required inputs. These changes prevent undocumented dataset-specific constants being interpreted as portable inference.
- Repaired the unavailable `windowed_hudson_fst` call using existing Hudson numerator/denominator and window-sum functions; documented the need to reconcile study values. Added a guard for non-estimable pairwise groups. This is an explicit compatibility repair, not a silent numerical-validation claim.
- Retained the corrected actual-week lag analysis alongside clearly distinguished sensitivity analyses. An index lag is not automatically a calendar-week lag.
- Added main separate-foreground selection-panel sources and retained direct-model convergence/supplementary sources. Omitted an alternate main-panel design not established as the V4 target.
- Kept conditional-rate-scaling qualifications and source-pair audits. Its provenance input hash keys now use workspace-relative paths; the source file hashes must remain unchanged.
- Excluded zero-byte scripts, redundant drafts, figures/workbooks, fitted models, restricted samples, genomic/read data, private helper modules and local audit notes from publication.

The source manifest distinguishes source-file hashes from released-file hashes. Syntax validation, synthetic guard tests and disclosure checks are reported separately from unavailable full-data numerical and visual validation.
