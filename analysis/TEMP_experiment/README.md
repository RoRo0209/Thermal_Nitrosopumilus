# Temperature incubation figure source

Set `THERMAL_PROJECT_ROOT` to a separate analysis workspace. Place `Day0_Day30_Day60_two_group_relative_abundance_organized.xlsx` in `TEMP_experiment/`. The input workbook is not included. The manuscript data-availability statement identifies the associated sequencing record; the prepared workbook and its complete raw-read processing provenance must be acquired separately.

Run `plot_temperature_incubation_main_figure.R` to write the V4 Fig. 3e,f components and supporting plots under `TEMP_experiment/Figure_Temperature_Incubation/`. Exports include `Figure_Temperature_Incubation_2panel_equal_width.pdf/svg/png`, `Supplementary_Fig_Temperature_Composition.*`, `Supplementary_Fig_Temperature_Drivers_heatmap.*`, and the corresponding source tables. Run `validate_temperature_figure_outputs.R` afterward to check table schemas, fraction closure, expected design dimensions and export presence.

Required workbook sheets and schema:

- `Genome Groups`: `Genome`, `Ecotype`, `Strain`, `Tree name`, `Included in two_group`.
- `Raw Data`: normalized sample identifier, day, temperature, replicate, genome, ecotype and relative abundance. The source reads the temperature and abundance columns positionally (columns 5 and 9); preserve the original schema and inspect the script's rename step before supplying a revised workbook.

The script selects HighT/LowT genomes, sums recruitment per sample and produces mean/SD summaries and composition plots. The day-zero sample is a baseline, not an independent replicate at each experimental temperature. Read-recruitment relative abundance is not an absolute growth-rate measure. No new hypothesis test is introduced by this port.

Requirements: R packages `readxl`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `scales`, `svglite`, `ragg`; optional `png` for validation image dimensions. The private `R_lib` injection was removed; install packages in the active R environment. Validation of this release was limited to source syntax; the workbook was not rerun and image appearance was not revalidated.
