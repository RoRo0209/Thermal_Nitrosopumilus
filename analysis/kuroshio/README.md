# Kuroshio regional validation: code without restricted data

**Kuroshio inputs are restricted and excluded from this repository.** Request the environmental and recruitment tables from their data owners. Set `THERMAL_PROJECT_ROOT` to a separate workspace with a `kuroshio/` directory. The scripts use that directory for all data and outputs.

The available sequence is:

1. `kuroshio_analysis_v3.R` joins `environmental_factor.xlsx` with `kuroshio_np_temp_abudance.xlsx`, then writes `output_tables/` and `output_figures/`. The environmental workbook was not located in the audited module; frozen downstream tables do not replace its provenance. Sheet 1 of the environmental file is read with one skipped row and 13 positional columns: `Station`, `Salinity`, `Temperature`, `Pressure`, `DO`, `Chla`, `Ammonium`, `Nitrate`, `Nitrite`, `Phosphate`, `Silicate`, `Density`, `CDOM`. The abundance workbook has a genome identifier column, an ecotype/group column and station/sample abundance columns. Station naming must match the parser in the script.
2. `Fig4B_ratio_vs_temperature_depth_filtered.R` reads `output_tables/meta_with_env_and_sum.csv` and optionally `output_tables/watermass_labels.csv`, applies the original depth filter and log-ratio pseudocount, and writes `Fig4B_ratio_vs_temperature_simplified.pdf` and `.png`. This is the available V4 Fig. 3c component despite the older filename.
3. `FigS6_spatial_gradient.R` reads the same source tables and writes `Supplementary_Fig_S6_depth_filtered.pdf` and `.png`. Its manuscript supplementary number may differ from the filename; see the release figure map.

The generated metadata table must retain the fields used by these plots: station/sample identity, `Temperature`, `HighT`, `LowT`, depth (`Depth_m` as derived by the script), and environmental/grouping fields. The code preserves the original temperature–salinity water-mass rules and assigns a fixed depth to names marked DCM. Verify those choices against the data-owner metadata rather than treating them as measured depths. Optional water-mass labels must have unique join keys.

The raw-read recruitment workflow producing the abundance workbook was not located here. Reference-genome identity, read filters, sample membership and the mapping from restricted raw files to workbook columns remain required external provenance.

Dependencies: R packages `readxl`, `tibble`, `dplyr`, `tidyr`, `stringr`, `forcats`, `purrr`, `vegan`, `ggplot2`, `ggrepel`, `cowplot`, `patchwork`, `scales`, `readr`. Source syntax checks only; no restricted-data rerun or numerical reproduction is claimed.
