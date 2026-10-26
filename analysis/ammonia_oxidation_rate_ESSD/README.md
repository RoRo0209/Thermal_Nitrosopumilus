# Published ocean nitrification database plots

Set `THERMAL_PROJECT_ROOT` to the separate analysis workspace and place `Database of nitrification and nitrifiers in the global ocean.xlsx` under `ammonia_oxidation_rate_ESSD/`. This workbook is a published source attributed to Tang and colleagues in the manuscript; obtain it from the publication's data record. It is not copied into this code release. The precise source URL/version/checksum must be reconciled with the manuscript data-availability record before a frozen reproduction run.

| Script | Role | Generated files |
|---|---|---|
| `plot_main_temperature_panels.py` | V4 Fig. 1b,c temperature panels | `Fig_main_panel_AOA_abundance_vs_temp.*`, `Fig_main_panel_AO_rate_vs_temp_by_NH4.*`, source-data CSVs and preview PNG |
| `plot_AOA_map_redrawn.py` | Supplementary sampling-map candidate | `Fig_AOA_map_redrawn.pdf`, `.png` |

The exact binned-abundance source for V4 Supplementary Fig. S2 is unresolved and is not replaced by a different temperature scatterplot here.

The main script reads sheets `ammonia oxidizer` and `volumetric ammonia oxidation`. Required source headers include `Depth (m)`, `Temperature (°C)`, `total AOA (copies L-1)`, `shallow clade AOA (copies L-1)`, `Rate (nmol N L-1 d-1)`, and the source ammonium/ammonia column. Map scripts also require `Latitude` and `Longitude`. The AOA pipeline fills missing total AOA abundance with shallow-clade abundance, retains positive abundance and depth at most 200 m, and log-transforms abundance. The rate script uses the source ammonium split and separate temperature associations. These calculations and filters are preserved from the available source; they do not establish an ecotype-specific causal effect on activity.

Run scripts by their absolute code paths. They switch to the workspace module and write there. Requirements: Python, NumPy, pandas, openpyxl, matplotlib, SciPy and PyMuPDF (`fitz`) for main-panel preview. The map needs Plotly and Kaleido with a compatible Chrome/Chromium renderer and geographic assets. Arial is preferred; a missing font may change layout. Final assembled manuscript figures are not rebuilt by these scripts alone.

Validation: Python syntax checked; no database rerun, numerical equality check or image-rendering validation is claimed.

## Conditional scaling and source-pair audits (S14 interpretation context)

Recovered additional code: `fixed_coefficient_scenario_20261016/calculate_and_plot.py`, `rate_link_revision_20261019/pair_audit/audit_pairs.py`, and `rate_link_revision_20261019/empirical_diagnostic/run_empirical_diagnostic.py`. Run by code path with `THERMAL_PROJECT_ROOT` set. They write into the corresponding workspace subdirectories.

The scaling script retains its original `exploratory_only_temperature_mismatch_unresolved_not_submission_ready` flag. It is a conditional fixed-coefficient exercise, not a validated future process-rate forecast. The later pair audit distinguishes matching metadata from proven shared samples/bottles. The diagnostic uses held-out study/grid structure without future extrapolation. These sources inform the manuscript's interpretation limits; their presence does not establish an exact final S14 redraw.

Inputs include the published observation workbook, the original Tolar event table, prepared model prediction grids and the intermediate `rate_abundance_overlay_preview_20261016` tables/summary, `coastal_rate_assessment_20261015/raw_selection_audit.tsv`, `high_temperature_search_20261016/ecotype_temperature_audit.json`, and `ecotype_function_20261002/paired_coastal_grids.tsv`. The original generation chain for all these intermediate files is not included. The scaling provenance manifest must use workspace-relative `input_hashes` keys while preserving the original file hashes. The two later audits additionally consume the scaling outputs. The final coastal-map rendering script depends on unavailable private figure-helper modules and was not packaged as a portable entry point.
