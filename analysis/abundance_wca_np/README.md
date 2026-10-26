# NP/WCA abundance and environmental models

These are available study scripts, ported to a separate analysis workspace. Set `THERMAL_PROJECT_ROOT` to an existing directory containing `abundance_wca_np/`. All relative inputs and outputs resolve there, not inside this code checkout. No observations, reads, reference sequences, metadata or fitted models are bundled. Numerical agreement with the V4 figures has not been rerun.

`calc_NP_WCA_ratio.R` is the source for the predicted map used in Fig. 1a. It reads `NP_WCA_final.xlsx` (first three columns: sample identifier, NP fraction, WCA fraction) and `genome_lon_lat.xlsx` (first three columns: sample identifier, longitude, latitude). The script joins these tables, normalizes fractions, joins gridded predictors, fits random forests and exports maps/tables under `NP_WCA_global_ratio/`, including `NP_WCA_global_map_RF_v2.pdf`. It also contains a future-grid projection section.

Required external input variables:

| Variable | Required file |
|---|---|
| `THERMAL_BASELINE_GRID` | `GLOBAL_grid_baseline_depthsurf_with_dist_km.tsv`; longitude/latitude, distance-to-coast and baseline predictor columns selected by the script |
| `THERMAL_BATHYMETRY_NC` | `terrain_characteristics.nc`; matched bathymetry grid |
| `THERMAL_FUTURE_GRID` | `GLOBAL_grid_future_SSP585_2080-2100_depthsurf.tsv`; predictor names and units matching the baseline |

These derived grids must be obtained or rebuilt from the documented Bio-ORACLE source. The grid-generation workflow and exact downloaded dataset versions are not supplied by this module. Do not interchange grids with different coordinates, units or variable definitions. The existing first-record duplicate removal and coordinate matching rules are retained; input join uniqueness should be audited for new datasets.

`supplementary_analysis_NP_WCA.R` uses the same workbook pair and baseline grid for partial correlations and GAMs. Its exact final supplementary-panel placement remains to be confirmed.

Available upstream classification scripts:

1. `run_oga_analysis_v3.sh <absolute-reference-FASTA>` reads `ABY89139.1_nucl_homologs.fna` and the tab-separated `abundance_matrix.csv` with two header rows, then writes `oga_results/`. Reference headers are `{ecotype}__{sequence_id}` for `HighT`, `LowT`, `NP_other`, `WCA`, `WCB`. The OGA files must be acquired separately. Both SRF and DCM records can be emitted; a surface-only dataset requires the appropriate downstream selection.
2. `coastal_amoA_diamond_v2.sh <absolute-reference-FASTA> <absolute-read-directory> [absolute-metadata-TSV]` expects paired `*_1.fastq.gz` and `*_2.fastq.gz`; optional metadata contains `sample`, `latitude`, `longitude`, `temperature`. It writes `coastal_amoA_results_v2/`. Restricted coastal data are not distributed. Existing output files are reused, so use a clean analysis workspace when changing the reference or inputs.

The conversion of classification outputs into the final model workbooks was not located as a complete, audited workflow. Reference-to-ecotype assignment provenance is an external requirement. This is a declared reproduction gap.

Dependencies: R packages `readxl`, `readr`, `dplyr`, `tidyr`, `ranger`, `ggplot2`, `scales`, `ncdf4`; supplementary analysis also `mgcv`, `ppcor`, `patchwork`. Upstream scripts need Bash, BLAST+, DIAMOND, Python with pandas/NumPy/matplotlib, and a nucleotide-to-protein translation command as used in the source. Linux/WSL is required for shell workflows. Use paths without whitespace for the legacy shell commands, which retain some unquoted arguments.

Validation: source syntax checks only; no study data run or model-performance validation is claimed.
