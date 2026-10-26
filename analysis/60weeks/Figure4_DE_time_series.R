#!/usr/bin/env Rscript
# Release path adapter: data and outputs live in a separate user workspace.
.release_root <- Sys.getenv("THERMAL_PROJECT_ROOT")
if (!nzchar(.release_root)) stop("Set THERMAL_PROJECT_ROOT to the analysis workspace")
.release_root <- normalizePath(.release_root, winslash = "/", mustWork = TRUE)
.release_input <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) stop("Set required input environment variable: ", name)
  normalizePath(value, winslash = "/", mustWork = TRUE)
}
setwd(file.path(.release_root, "60weeks"))

# ============================================================
# Figure 4 Panel D & E: Time series (Temperature + Ecotype ratio)
# Modified from 60weeks_pub.R Panel A
# Split dual y-axis into two independent panels
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ---- 1. Paths ----
# Use relative paths from current working directory
indir   <- "out_60weeks_by_station_mgonly"
in_csv  <- file.path(indir, "60weeks_station_week_table_mgonly.csv")

# ---- 2. Read & prepare data ----
dat <- read_csv(in_csv, show_col_types = FALSE) %>%
  mutate(
    station     = factor(station, levels = c("E", "W", "Z")),
    week_num    = as.integer(week_num),
    Temperature = as.numeric(Temperature),
    HighT       = as.numeric(HighT),
    LowT        = as.numeric(LowT)
  ) %>%
  filter(is.finite(Temperature), is.finite(HighT), is.finite(LowT))

# Calculate ratio if not present or recalculate to ensure consistency
if (!"ratio" %in% colnames(dat) || any(is.na(dat$ratio))) {
  dat <- dat %>%
    mutate(ratio = log10((HighT + 1e-6) / (LowT + 1e-6)))
} else {
  dat <- dat %>%
    mutate(ratio = as.numeric(ratio))
}

# Print data ranges
cat("\n=== Data Summary ===\n")
cat("Temperature range:", min(dat$Temperature, na.rm = TRUE), "to",
    max(dat$Temperature, na.rm = TRUE), "°C\n")
cat("Ratio range:", min(dat$ratio, na.rm = TRUE), "to",
    max(dat$ratio, na.rm = TRUE), "\n\n")

# ---- 3. Style settings ----
# ISME-style colorblind-safe colors
temp_color <- "#2C5F8D"           # Deep blue for temperature
temp_ribbon <- "#B8D4E8"          # Light blue for temperature range
ratio_mean_color <- "#C41E3A"     # Deep red for mean ratio
ratio_ribbon <- "#F4C2C2"         # Light red/pink for ratio range
ref_line_color <- "#666666"       # Dark grey for reference line

# Base theme following ISME/Nature style
base_theme <- theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.7),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "none",  # Remove legend for cleaner look
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    plot.margin = margin(6, 10, 6, 6)
  )

# ---- 4. Prepare summary data ----
# For Panel D: mean temperature with min-max envelope
dat_temp_summary <- dat %>%
  group_by(week_num) %>%
  summarize(
    mean_temp = mean(Temperature, na.rm = TRUE),
    min_temp  = min(Temperature, na.rm = TRUE),
    max_temp  = max(Temperature, na.rm = TRUE),
    .groups = "drop"
  )

# For Panel E: mean ratio with station range
dat_ratio_summary <- dat %>%
  group_by(week_num) %>%
  summarize(
    mean_ratio = mean(ratio, na.rm = TRUE),
    min_ratio  = min(ratio, na.rm = TRUE),
    max_ratio  = max(ratio, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 5. Panel D: Seasonal temperature envelope ----
pD <- ggplot() +
  # Station range ribbon (light blue shadow)
  geom_ribbon(data = dat_temp_summary,
              aes(x = week_num, ymin = min_temp, ymax = max_temp),
              fill = temp_ribbon, alpha = 0.5) +
  # Mean temperature line (deep blue)
  geom_line(data = dat_temp_summary,
            aes(x = week_num, y = mean_temp),
            color = temp_color, linewidth = 1.2) +
  geom_point(data = dat_temp_summary,
             aes(x = week_num, y = mean_temp),
             color = temp_color, size = 2.2, shape = 19) +
  scale_x_continuous(breaks = seq(0, 60, 10), expand = c(0.02, 0)) +
  scale_y_continuous(breaks = seq(5, 30, 5)) +
  labs(
    title = "D  Seasonal temperature cycles",
    x = NULL,
    y = "Temperature (°C)"
  ) +
  base_theme +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# ---- 6. Panel E: Ecotype ratio with station range shadow ----
pE <- ggplot() +
  # y=0 reference line
  geom_hline(yintercept = 0, linetype = 2, color = ref_line_color, linewidth = 0.8) +
  # Station range ribbon (light red/pink shadow) - elegant and clean!
  geom_ribbon(data = dat_ratio_summary,
              aes(x = week_num, ymin = min_ratio, ymax = max_ratio),
              fill = ratio_ribbon, alpha = 0.5) +
  # Mean ratio line (deep red, emphasized)
  geom_line(data = dat_ratio_summary,
            aes(x = week_num, y = mean_ratio),
            color = ratio_mean_color, linewidth = 1.5) +
  geom_point(data = dat_ratio_summary,
             aes(x = week_num, y = mean_ratio),
             color = ratio_mean_color, size = 3, shape = 19) +
  scale_x_continuous(breaks = seq(0, 60, 10), expand = c(0.02, 0)) +
  scale_y_continuous(breaks = seq(-3.5, 0, 0.5)) +
  labs(
    title = "E  Persistent LowT dominance across seasonal cycle",
    x = "Week",
    y = expression(log[10](HighT/LowT))
  ) +
  base_theme

# ---- 7. Combine panels ----
fig <- pD / pE +
  plot_layout(heights = c(0.75, 1.25))

# ---- 8. Save outputs ----
out_pdf <- file.path(indir, "Fig4_DE_ISME_style.pdf")
out_png <- file.path(indir, "Fig4_DE_ISME_style.png")

# Adjusted dimensions: wide aspect ratio (2:1) for horizontal layout
ggsave(out_pdf, fig, width = 8, height = 4, units = "in")
ggsave(out_png, fig, width = 8, height = 4, units = "in", dpi = 600)  # ISME standard: 600 dpi

cat("✓ Figures saved (ISME style):\n")
cat("  -", out_pdf, "\n")
cat("  -", out_png, "(8 × 4 in, 600 dpi)\n")
cat("\nDone!\n")
