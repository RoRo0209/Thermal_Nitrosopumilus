# Release path adapter: data and outputs live in a separate user workspace.
.release_root <- Sys.getenv("THERMAL_PROJECT_ROOT")
if (!nzchar(.release_root)) stop("Set THERMAL_PROJECT_ROOT to the analysis workspace")
.release_root <- normalizePath(.release_root, winslash = "/", mustWork = TRUE)
.release_input <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) stop("Set required input environment variable: ", name)
  normalizePath(value, winslash = "/", mustWork = TRUE)
}
setwd(.release_root)

options(stringsAsFactors = FALSE)

# Install dependencies in the active R environment; no private library injection.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
})

input_file <- "TEMP_experiment/Day0_Day30_Day60_two_group_relative_abundance_organized.xlsx"
output_dir <- "TEMP_experiment/Figure_Temperature_Incubation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

high_col <- "#D75A4A"
low_col <- "#3E78A8"
neutral_col <- "#333333"
grid_col <- "#D9D9D9"
font_family <- "Arial"

theme_paper <- function(base_size = 6.7) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.title = element_text(size = 6.7, colour = "black"),
      axis.text = element_text(size = 6.1, colour = "black"),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.32, colour = "black"),
      axis.ticks.length = unit(1.4, "mm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 6.1),
      legend.key.height = unit(3.0, "mm"),
      legend.key.width = unit(4.2, "mm"),
      plot.margin = margin(5, 5, 5, 5, unit = "pt"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

groups <- read_excel(input_file, sheet = "Genome Groups") %>%
  transmute(
    Genome = as.character(Genome),
    Metadata_ecotype = as.character(Ecotype),
    Strain = as.character(Strain),
    Tree_name = as.character(`Tree name`),
    Included = as.character(`Included in two_group`)
  ) %>%
  filter(Metadata_ecotype %in% c("HighT", "LowT")) %>%
  distinct(Genome, .keep_all = TRUE)

raw_input <- read_excel(input_file, sheet = "Raw Data")
raw <- raw_input %>%
  transmute(
    Sample = as.character(`Normalized sample`),
    Day = as.numeric(Day),
    Temperature = as.numeric(.data[[names(raw_input)[5]]]),
    Replicate = as.numeric(Replicate),
    Genome = as.character(Genome),
    Workbook_ecotype = as.character(Ecotype),
    Relative_abundance = as.numeric(.data[[names(raw_input)[9]]])
  ) %>%
  filter(Genome %in% groups$Genome) %>%
  left_join(groups, by = "Genome")

if (any(is.na(raw$Metadata_ecotype))) stop("Some genomes did not match genome metadata.")
if (any(raw$Workbook_ecotype != raw$Metadata_ecotype, na.rm = TRUE)) stop("Ecotype mismatch between Raw Data and Genome Groups.")

sample_levels <- c(
  "T0",
  paste0("T12_D30_R", 1:3), paste0("T24_D30_R", 1:3), paste0("T30_D30_R", 1:3),
  paste0("T12_D60_R", 1:3), paste0("T24_D60_R", 1:3), paste0("T30_D60_R", 1:3)
)

summary_wide <- raw %>%
  group_by(Sample, Day, Temperature, Replicate, Metadata_ecotype) %>%
  summarise(Abundance = sum(Relative_abundance, na.rm = TRUE), .groups = "drop") %>%
  complete(Sample, Metadata_ecotype = c("HighT", "LowT"), fill = list(Abundance = 0)) %>%
  pivot_wider(names_from = Metadata_ecotype, values_from = Abundance, values_fill = 0) %>%
  mutate(
    Time = paste0("Day", Day),
    Total_two_group = HighT + LowT,
    HighT_fraction = if_else(Total_two_group > 0, HighT / Total_two_group, NA_real_),
    LowT_fraction = if_else(Total_two_group > 0, LowT / Total_two_group, NA_real_),
    Sample_order = match(Sample, sample_levels)
  ) %>%
  arrange(Sample_order)

if (!identical(summary_wide$Sample, sample_levels)) stop("Unexpected sample names or order.")
if (nrow(summary_wide) != 19) stop("Expected 19 samples.")

summary_export <- summary_wide %>%
  transmute(
    Sample,
    Temperature,
    Time,
    Replicate,
    HighT_abundance = HighT,
    LowT_abundance = LowT,
    Total_two_group,
    HighT_fraction,
    LowT_fraction
  )

write.table(
  summary_export,
  file.path(output_dir, "SourceData_Temperature_Incubation_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

group_stats <- summary_export %>%
  filter(Time != "Day0") %>%
  group_by(Time, Temperature) %>%
  summarise(
    n = n(),
    HighT_fraction_mean = mean(HighT_fraction),
    HighT_fraction_SD = sd(HighT_fraction),
    Total_two_group_mean = mean(Total_two_group),
    Total_two_group_SD = sd(Total_two_group),
    .groups = "drop"
  )
write.table(group_stats, file.path(output_dir, "GroupStatistics_Temperature_Incubation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Shared numeric x positions preserve replicate bars while leaving visible gaps.
sample_positions <- tibble(
  Sample = sample_levels,
  xpos = c(
    1,
    2.90, 3.25, 3.60, 4.70, 5.05, 5.40, 6.50, 6.85, 7.20,
    9.10, 9.45, 9.80, 10.90, 11.25, 11.60, 12.70, 13.05, 13.40
  )
)
group_centres <- c(1, 3.25, 5.05, 6.85, 9.45, 11.25, 13.05)
group_labels <- c("Day 0", "12", "24", "30", "12", "24", "30")

plot_df <- summary_wide %>% left_join(sample_positions, by = "Sample")

# Panel A: composition within target ecotypes.
panel_a_df <- plot_df %>%
  select(Sample, xpos, HighT_fraction, LowT_fraction) %>%
  pivot_longer(c(HighT_fraction, LowT_fraction), names_to = "Ecotype", values_to = "Fraction") %>%
  mutate(
    Ecotype = recode(Ecotype, HighT_fraction = "HighT", LowT_fraction = "LowT"),
    Ecotype = factor(Ecotype, levels = c("LowT", "HighT"))
  )

p_a <- ggplot(panel_a_df, aes(x = xpos, y = Fraction * 100, fill = Ecotype)) +
  geom_col(width = 0.27, colour = "white", linewidth = 0.18) +
  geom_vline(xintercept = c(2.05, 8.15), colour = grid_col, linewidth = 0.38) +
  annotate("text", x = 1.0, y = 105.5, label = "Day 0", family = font_family, size = 2.25) +
  annotate("text", x = 5.05, y = 105.5, label = "Day 30", family = font_family, size = 2.25) +
  annotate("text", x = 11.25, y = 105.5, label = "Day 60", family = font_family, size = 2.25) +
  scale_fill_manual(values = c(LowT = low_col, HighT = high_col), breaks = c("HighT", "LowT")) +
  scale_x_continuous(breaks = group_centres, labels = group_labels, expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 109), expand = c(0, 0)) +
  labs(x = expression(Temperature~(degree*C)), y = "Relative fraction within\ntarget ecotypes (%)") +
  coord_cartesian(clip = "off") +
  theme_paper() +
  theme(
    legend.position = "top", legend.direction = "horizontal", legend.background = element_blank(),
    legend.margin = margin(0, 0, -2, 0), legend.box.spacing = unit(0, "pt"),
    axis.title.x = element_text(margin = margin(t = 3)),
    plot.margin = margin(8, 5, 5, 5, unit = "pt")
  )

# Main figure: time trajectories. Temperature is encoded consistently across panels.
temp_cols <- c(`12` = "#3B6FB6", `24` = "#D28E00", `30` = "#C74440")
baseline_fraction <- summary_export$HighT_fraction[summary_export$Time == "Day0"] * 100
baseline_total <- summary_export$Total_two_group[summary_export$Time == "Day0"]

trajectory_raw <- summary_export %>%
  filter(Time %in% c("Day30", "Day60")) %>%
  mutate(
    Temperature_f = factor(as.character(Temperature), levels = c("12", "24", "30")),
    x_jitter = if_else(Time == "Day30", 30, 60) + (Replicate - 2) * 1.05,
    HighT_pct = HighT_fraction * 100
  )

trajectory_stats <- trajectory_raw %>%
  group_by(Temperature_f, Temperature, Time) %>%
  summarise(
    Day = if_else(first(Time) == "Day30", 30, 60),
    HighT_mean = mean(HighT_pct), HighT_sd = sd(HighT_pct),
    Total_mean = mean(Total_two_group), Total_sd = sd(Total_two_group),
    .groups = "drop"
  ) %>%
  select(-Time) %>%
  bind_rows(tibble(
    Temperature_f = factor(c("12", "24", "30"), levels = c("12", "24", "30")),
    Temperature = c(12, 24, 30), Day = 0,
    HighT_mean = baseline_fraction, HighT_sd = NA_real_,
    Total_mean = baseline_total, Total_sd = NA_real_
  )) %>%
  arrange(Temperature, Day)

end_labels <- trajectory_stats %>%
  filter(Day == 60) %>%
  mutate(Label = paste0(Temperature, "*degree*C"))

# HighT-composition trajectories through time (displayed as Panel B).
p_traj_a <- ggplot() +
  geom_hline(yintercept = 50, colour = grid_col, linewidth = 0.40, linetype = "dotted") +
  geom_line(data = trajectory_stats,
            aes(x = Day, y = HighT_mean, colour = Temperature_f, group = Temperature_f),
            linewidth = 0.78, lineend = "round") +
  geom_errorbar(data = filter(trajectory_stats, Day > 0),
                aes(x = Day, ymin = pmax(0, HighT_mean - HighT_sd),
                    ymax = pmin(100, HighT_mean + HighT_sd), colour = Temperature_f),
                width = 1.6, linewidth = 0.48) +
  geom_point(data = trajectory_raw,
             aes(x = x_jitter, y = HighT_pct, colour = Temperature_f),
             size = 1.45, alpha = 0.82) +
  geom_point(data = filter(trajectory_stats, Day > 0),
             aes(x = Day, y = HighT_mean, colour = Temperature_f),
             shape = 21, fill = "white", size = 1.75, stroke = 0.72) +
  geom_point(aes(x = 0, y = baseline_fraction), shape = 21, fill = "white",
             colour = neutral_col, size = 1.85, stroke = 0.72) +
  geom_text(data = end_labels,
            aes(x = 61.0, y = HighT_mean, label = Label, colour = Temperature_f),
            parse = TRUE, hjust = 0, family = font_family, fontface = "bold", size = 1.75) +
  annotate("text", x = 39, y = 23.5, label = "LowT-dominated",
           family = font_family, size = 1.45, colour = temp_cols[["12"]]) +
  annotate("text", x = 40, y = 70.5, label = "Transient enrichment\nthen coexistence",
           family = font_family, size = 1.45, colour = temp_cols[["24"]]) +
  annotate("text", x = 42, y = 93.0, label = "HighT-dominated",
           family = font_family, size = 1.45, colour = temp_cols[["30"]]) +
  scale_colour_manual(values = temp_cols) +
  scale_x_continuous(breaks = c(0, 30, 60), labels = c("Day 0", "Day 30", "Day 60"),
                     limits = c(-2, 69), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 103), expand = c(0, 0)) +
  labs(x = NULL, y = "HighT fraction within target ecotypes (%)") +
  coord_cartesian(clip = "off") +
  theme_paper(base_size = 7.0) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 5.2, margin = margin(r = 3)),
    axis.text = element_text(size = 5.4),
    plot.margin = margin(13, 18, 13, 5, unit = "pt")
  )

# Total target-ecotype abundance trajectories (displayed as Panel A).
log_guides <- 10^seq(-4, 0, 1)
p_traj_b <- ggplot() +
  geom_hline(yintercept = log_guides, colour = "#ECECEC", linewidth = 0.32) +
  geom_line(data = trajectory_stats,
            aes(x = Day, y = Total_mean, colour = Temperature_f, group = Temperature_f),
            linewidth = 0.72, lineend = "round") +
  geom_errorbar(data = filter(trajectory_stats, Day > 0),
                aes(x = Day, ymin = pmax(Total_mean - Total_sd, 7e-5),
                    ymax = Total_mean + Total_sd, colour = Temperature_f),
                width = 1.8, linewidth = 0.45) +
  geom_point(data = trajectory_raw,
             aes(x = x_jitter, y = Total_two_group, colour = Temperature_f),
             size = 1.35, alpha = 0.82) +
  geom_point(data = filter(trajectory_stats, Day > 0),
             aes(x = Day, y = Total_mean, colour = Temperature_f),
             shape = 21, fill = "white", size = 1.65, stroke = 0.68) +
  geom_point(aes(x = 0, y = baseline_total), shape = 21, fill = "white",
             colour = neutral_col, size = 1.75, stroke = 0.70) +
  scale_colour_manual(values = temp_cols) +
  scale_x_continuous(breaks = c(0, 30, 60), labels = c("Day 0", "Day 30", "Day 60"),
                     limits = c(-2, 63), expand = c(0, 0)) +
  scale_y_log10(
    breaks = log_guides,
    labels = function(x) parse(text = sprintf("10^{%d}", round(log10(x)))),
    limits = c(7e-5, 3.3), expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Relative abundance of target\nNitrosopumilus ecotypes (%)") +
  theme_paper(base_size = 7.0) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 5.2, margin = margin(r = 3)),
    axis.text = element_text(size = 5.4),
    plot.margin = margin(13, 5, 13, 5, unit = "pt")
  )

main_figure <- (p_traj_b | p_traj_a) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = font_family, face = "bold", size = 8.2, colour = "black"),
        plot.tag.position = c(0.22, 1.06))

set_patchwork_panel_size <- function(plot, width_mm, height_mm) {
  capture <- agg_capture(width = 183, height = 78, units = "mm", res = 96)
  on.exit(dev.off(), add = TRUE)
  grob <- patchworkGrob(plot)
  panel_layout <- grepl("^panel-[0-9]+$", grob$layout$name)
  panel_columns <- unique(grob$layout$l[panel_layout])
  panel_rows <- unique(grob$layout$t[panel_layout])
  grob$widths[panel_columns] <- grid::unit(width_mm, "mm")
  grob$heights[panel_rows] <- grid::unit(height_mm, "mm")
  list(
    grob = grob,
    width_in = grid::convertWidth(sum(grob$widths), "in", valueOnly = TRUE),
    height_in = grid::convertHeight(sum(grob$heights), "in", valueOnly = TRUE)
  )
}

main_sized <- set_patchwork_panel_size(main_figure, width_mm = 29.425, height_mm = 33.264)
main_width <- main_sized$width_in
main_height <- main_sized$height_in
ggsave(file.path(output_dir, "Figure_Temperature_Incubation_2panel_equal_width.pdf"), main_sized$grob,
       device = cairo_pdf, width = main_width, height = main_height, units = "in", family = font_family)
ggsave(file.path(output_dir, "Figure_Temperature_Incubation_2panel_equal_width.svg"), main_sized$grob,
       device = svglite, width = main_width, height = main_height, units = "in", bg = "white")
ggsave(file.path(output_dir, "Figure_Temperature_Incubation_2panel_equal_width.png"), main_sized$grob,
       device = agg_png, width = main_width, height = main_height, units = "in", dpi = 600, background = "white")

# The former stacked-composition panel is retained as a standalone supplementary figure.
supp_composition <- p_a +
  labs(x = expression(Temperature~(degree*C)), y = "Relative fraction within target ecotypes (%)") +
  theme(plot.margin = margin(5, 7, 5, 7, unit = "pt"))
supp_comp_width <- 150 / 25.4
supp_comp_height <- 70 / 25.4
ggsave(file.path(output_dir, "Supplementary_Fig_Temperature_Composition.pdf"), supp_composition,
       device = cairo_pdf, width = supp_comp_width, height = supp_comp_height, units = "in", family = font_family)
ggsave(file.path(output_dir, "Supplementary_Fig_Temperature_Composition.svg"), supp_composition,
       device = svglite, width = supp_comp_width, height = supp_comp_height, units = "in", bg = "white")
ggsave(file.path(output_dir, "Supplementary_Fig_Temperature_Composition.png"), supp_composition,
       device = agg_png, width = supp_comp_width, height = supp_comp_height, units = "in", dpi = 600, background = "white")

# Supplementary heatmap: reproducible selection of abundant treatment drivers.
treatment_means <- raw %>%
  filter(Day %in% c(30, 60)) %>%
  group_by(Genome, Metadata_ecotype, Strain, Day, Temperature) %>%
  summarise(Treatment_mean = mean(Relative_abundance), .groups = "drop")

driver_scores <- treatment_means %>%
  group_by(Genome, Metadata_ecotype, Strain) %>%
  summarise(Max_treatment_mean = max(Treatment_mean), .groups = "drop") %>%
  left_join(
    raw %>% group_by(Genome) %>% summarise(Detected_samples = sum(Relative_abundance > 0), Mean_all_samples = mean(Relative_abundance), .groups = "drop"),
    by = "Genome"
  ) %>%
  filter(Detected_samples >= 3) %>%
  group_by(Metadata_ecotype) %>%
  arrange(desc(Max_treatment_mean), desc(Mean_all_samples), .by_group = TRUE) %>%
  slice_head(n = 6) %>%
  ungroup()

if (!all(c("HighT", "LowT") %in% driver_scores$Metadata_ecotype)) stop("Driver selection lost an ecotype.")

nonzero <- raw$Relative_abundance[raw$Relative_abundance > 0]
pseudocount <- min(nonzero) / 2

selected_drivers <- driver_scores %>%
  mutate(Pseudocount = pseudocount) %>%
  arrange(factor(Metadata_ecotype, levels = c("HighT", "LowT")), desc(Max_treatment_mean))
write.table(selected_drivers, file.path(output_dir, "Selected_driver_genomes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

heat_df <- raw %>%
  filter(Genome %in% selected_drivers$Genome) %>%
  mutate(
    Sample = factor(Sample, levels = sample_levels),
    Ecotype = factor(Metadata_ecotype, levels = c("HighT", "LowT")),
    Genome_accession = sub("^((GCA|GCF)_[0-9]+\\.[0-9]+).*$", "\\1", Genome),
    Genome_label = if_else(!is.na(Strain) & Strain != "" & Strain != "NA", Strain, Genome_accession),
    Genome_label = if_else(nchar(Genome_label) > 22, paste0(substr(Genome_label, 1, 21), "..."), Genome_label),
    log10_abundance = log10(Relative_abundance + pseudocount)
  )

row_order <- selected_drivers %>%
  mutate(Genome_accession = sub("^((GCA|GCF)_[0-9]+\\.[0-9]+).*$", "\\1", Genome),
         Genome_label = if_else(!is.na(Strain) & Strain != "" & Strain != "NA", Strain, Genome_accession),
         Genome_label = if_else(nchar(Genome_label) > 22, paste0(substr(Genome_label, 1, 21), "..."), Genome_label)) %>%
  arrange(factor(Metadata_ecotype, levels = c("LowT", "HighT")), Max_treatment_mean) %>%
  pull(Genome_label)
heat_df$Genome_label <- factor(heat_df$Genome_label, levels = unique(row_order))

write.table(
  heat_df %>% transmute(Sample = as.character(Sample), Temperature, Time = paste0("Day", Day), Replicate,
                        Genome, Strain, Ecotype = as.character(Ecotype), Relative_abundance, log10_abundance, Pseudocount = pseudocount),
  file.path(output_dir, "SourceData_Temperature_Drivers_heatmap.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)

heat_min <- floor(min(heat_df$log10_abundance))
heat_max <- ceiling(max(heat_df$log10_abundance))

p_heat <- ggplot(heat_df, aes(x = Sample, y = Genome_label, fill = log10_abundance)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  facet_grid(Ecotype ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_vline(xintercept = c(1.5, 10.5), colour = "white", linewidth = 1.0) +
  scale_fill_gradientn(
    colours = c("#172A3A", "#315A7D", "#69A2B8", "#D7E7D0", "#F2C14E"),
    limits = c(heat_min, heat_max), oob = squish,
    name = expression(log[10](relative~abundance + pseudocount))
  ) +
  scale_x_discrete(labels = c(
    "T0" = "D0",
    setNames(paste0("12-R", 1:3), paste0("T12_D30_R", 1:3)),
    setNames(paste0("24-R", 1:3), paste0("T24_D30_R", 1:3)),
    setNames(paste0("30-R", 1:3), paste0("T30_D30_R", 1:3)),
    setNames(paste0("12-R", 1:3), paste0("T12_D60_R", 1:3)),
    setNames(paste0("24-R", 1:3), paste0("T24_D60_R", 1:3)),
    setNames(paste0("30-R", 1:3), paste0("T30_D60_R", 1:3))
  )) +
  labs(x = "Sample (temperature and replicate)", y = NULL,
       title = "Genomes driving the temperature response",
       subtitle = "Day 0                     Day 30                                      Day 60") +
  theme_minimal(base_size = 7, base_family = font_family) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.1, colour = "black"),
    axis.text.y = element_text(size = 6.4, colour = "black"),
    axis.title.x = element_text(size = 7, margin = margin(t = 5)),
    strip.placement = "outside",
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 7, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6.2),
    plot.title = element_text(size = 8.2, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 6.4, colour = "#444444"),
    plot.margin = margin(6, 7, 6, 6, unit = "pt"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

heat_width <- 160 / 25.4
heat_height <- 90 / 25.4
ggsave(file.path(output_dir, "Supplementary_Fig_Temperature_Drivers_heatmap.pdf"), p_heat,
       device = cairo_pdf, width = heat_width, height = heat_height, units = "in", family = font_family)
ggsave(file.path(output_dir, "Supplementary_Fig_Temperature_Drivers_heatmap.svg"), p_heat,
       device = svglite, width = heat_width, height = heat_height, units = "in", bg = "white")
ggsave(file.path(output_dir, "Supplementary_Fig_Temperature_Drivers_heatmap.png"), p_heat,
       device = agg_png, width = heat_width, height = heat_height, units = "in", dpi = 600, background = "white")

legend_text <- c(
  "# Figure legend",
  "",
  "**Temperature-dependent temporal trajectories of Nitrosopumilus HighT and LowT ecotypes during a 60-day coastal seawater incubation.** ",
  "**A,** Combined metagenomic relative abundance of HighT and LowT target ecotypes through time at 12, 24 and 30 °C. All trajectories share the single Day 0 sample. Small symbols show independent biological replicates on Days 30 and 60 (n = 3 per temperature and time point); large open symbols and error bars show mean ± s.d. Lines connect treatment means only, and the y axis is logarithmic. Relative abundance is compositional and is not an absolute abundance or activity measurement. ",
  "**B,** HighT fraction within the mapped target ecotypes through time. Symbols and error bars are defined as in A. The dotted horizontal line denotes equal HighT and LowT representation.",
  "",
  "# Supplementary composition legend",
  "",
  "**Relative composition of HighT and LowT genomes within target Nitrosopumilus ecotypes.** Day 0 is represented by one sample; Days 30 and 60 show three independent biological replicates for each temperature. Each bar sums to 100%.",
  "",
  "# Supplementary heatmap legend",
  "",
  "**Representative genomes driving the temperature response.** Six genomes per ecotype were selected by the highest treatment-level mean relative abundance after requiring detection in at least three samples. Colours show log10(relative abundance + pseudocount), where the pseudocount is one-half of the smallest non-zero genome abundance in the full two-ecotype dataset."
)
writeLines(legend_text, file.path(output_dir, "Figure_Temperature_Incubation_legend.md"), useBytes = TRUE)

qa_text <- c(
  "# Figure QA notes",
  "",
  paste0("- Input workbook: ", basename(input_file)),
  paste0("- Samples: ", nrow(summary_export), " (Day 0 n=1; each Day × temperature treatment n=3)"),
  paste0("- Target genomes: ", n_distinct(raw$Genome), " (HighT ", n_distinct(raw$Genome[raw$Metadata_ecotype == "HighT"]), "; LowT ", n_distinct(raw$Genome[raw$Metadata_ecotype == "LowT"]), ")"),
  "- Aggregation: genome relative abundances summed within metadata-defined ecotypes for each sample.",
  "- Fractions: HighT or LowT divided by HighT + LowT; all values finite and bars sum to 100%.",
  "- Statistics: mean ± s.d.; no inferential tests were added because none were requested.",
  "- Day 0 is a visual baseline only and is excluded from temperature statistics.",
  "- Panel A uses log10 scaling; all target totals are positive.",
  "- Main-panel layout: total abundance at left and HighT fraction at right, with equal panel widths.",
  "- The former stacked-composition panel is retained as a standalone supplementary figure.",
  paste0("- Heatmap pseudocount: ", format(pseudocount, scientific = TRUE, digits = 6)),
  "- Interpretation boundary: CoverM relative abundance is compositional and does not establish absolute growth, nitrification rate or 15N activity.",
  "- Exports: editable PDF/SVG and 600-dpi PNG; white background; Arial typeface requested."
)
writeLines(qa_text, file.path(output_dir, "Figure_QA_notes.md"), useBytes = TRUE)

cat("Created outputs in", normalizePath(output_dir, winslash = "/"), "\n")
cat("Day0 HighT baseline (%):", sprintf("%.4f", baseline_fraction), "\n")
cat("Heatmap pseudocount:", format(pseudocount, scientific = TRUE), "\n")
print(selected_drivers %>% select(Metadata_ecotype, Genome, Strain, Max_treatment_mean, Detected_samples))
