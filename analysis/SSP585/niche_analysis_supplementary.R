# Public code port: data are supplied separately; see the module README.
THERMAL_PROJECT_ROOT <- Sys.getenv("THERMAL_PROJECT_ROOT")
if (!nzchar(THERMAL_PROJECT_ROOT) || !dir.exists(THERMAL_PROJECT_ROOT)) {
  stop("Set THERMAL_PROJECT_ROOT to an existing, separate analysis workspace.", call. = FALSE)
}
THERMAL_PROJECT_ROOT <- normalizePath(THERMAL_PROJECT_ROOT, winslash = "/", mustWork = TRUE)
thermal_path <- function(...) file.path(THERMAL_PROJECT_ROOT, ...)
thermal_required_env <- function(name, directory = FALSE) {
  value <- Sys.getenv(name)
  if (!nzchar(value) || (directory && !dir.exists(value)) || (!directory && !file.exists(value))) {
    stop("Missing required input: set ", name, " to an existing ", if (directory) "directory." else "file.", call. = FALSE)
  }
  normalizePath(value, winslash = "/", mustWork = TRUE)
}
setwd(thermal_path("SSP585"))


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- thermal_path("SSP585")
out_dir <- file.path(base_dir, "niche_analysis_output")
lowt_file <- file.path(base_dir, "NP_lowT/AOA_points_with_BioORACLE_nearshore_filtered.tsv")
hight_file <- file.path(base_dir, "NP_highT/AOA_points_with_BioORACLE_nearshore_filtered.tsv")

env_vars <- c("temp", "sal", "chl", "o2", "no3", "po4", "si")
labels <- c(
  temp = "Temperature", sal = "Salinity", chl = "Chlorophyll-a",
  o2 = "Oxygen", no3 = "Nitrate", po4 = "Phosphate", si = "Silicate",
  Latitude = "Latitude"
)
col_lowt <- "#8EC9D8"
col_hight <- "#E89A72"

lowt <- read_tsv(lowt_file, show_col_types = FALSE) |>
  mutate(ecotype = "LowT") |>
  filter(if_all(all_of(env_vars), is.finite))
hight <- read_tsv(hight_file, show_col_types = FALSE) |>
  mutate(ecotype = "HighT") |>
  filter(if_all(all_of(env_vars), is.finite))
dat <- bind_rows(lowt, hight)

weighted_sd <- function(x, w) {
  mu <- weighted.mean(x, w)
  sqrt(sum(w * (x - mu)^2) / sum(w))
}

mu_lowt <- weighted.mean(lowt$temp, lowt$AOA_ab)
mu_hight <- weighted.mean(hight$temp, hight$AOA_ab)
sd_lowt <- weighted_sd(lowt$temp, lowt$AOA_ab)
sd_hight <- weighted_sd(hight$temp, hight$AOA_ab)

dens_lowt <- density(
  lowt$temp, weights = lowt$AOA_ab / sum(lowt$AOA_ab),
  from = -2, to = 32, n = 512, adjust = 1.5
)
dens_hight <- density(
  hight$temp, weights = hight$AOA_ab / sum(hight$AOA_ab),
  from = -2, to = 32, n = 512, adjust = 1.5
)
dens_df <- bind_rows(
  data.frame(temp = dens_lowt$x, density = dens_lowt$y, ecotype = "LowT"),
  data.frame(temp = dens_hight$x, density = dens_hight$y, ecotype = "HighT")
)

base_theme <- theme_classic(base_size = 8.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(5, 7, 5, 5)
  )

p_a <- ggplot(dens_df, aes(temp, density, fill = ecotype, colour = ecotype)) +
  geom_area(alpha = 0.32, position = "identity") +
  geom_line(linewidth = 0.75) +
  geom_vline(xintercept = c(mu_lowt, mu_hight),
             colour = c(col_lowt, col_hight), linetype = "dashed", linewidth = 0.5) +
  annotate(
    "text", x = 1, y = max(dens_df$density) * 0.96, hjust = 0, vjust = 1,
    size = 2.3,
    label = sprintf(
      "LowT: %.1f ± %.1f °C\nHighT: %.1f ± %.1f °C\nSchoener's D = 0.272",
      mu_lowt, sd_lowt, mu_hight, sd_hight
    )
  ) +
  scale_fill_manual(values = c(LowT = col_lowt, HighT = col_hight)) +
  scale_colour_manual(values = c(LowT = col_lowt, HighT = col_hight)) +
  coord_cartesian(xlim = c(0, 32), expand = FALSE) +
  labs(x = "Sea-surface temperature (°C)", y = "Abundance-weighted density") +
  base_theme + theme(legend.position = "none")

levins_results <- list()
for (v in c(env_vars, "Latitude")) {
  all_vals <- if (v == "Latitude") dat$Latitude else dat[[v]]
  breaks <- seq(min(all_vals), max(all_vals), length.out = 21)
  for (eco in c("LowT", "HighT")) {
    sub <- filter(dat, ecotype == eco)
    vals <- if (v == "Latitude") sub$Latitude else sub[[v]]
    bins <- cut(vals, breaks, include.lowest = TRUE)
    p <- tapply(sub$AOA_ab, bins, sum, default = 0)
    p[is.na(p)] <- 0
    p <- p / sum(p)
    b <- 1 / sum(p[p > 0]^2)
    levins_results[[length(levins_results) + 1]] <- data.frame(
      var = v, ecotype = eco, B_std = (b - 1) / 19
    )
  }
}
levins_df <- bind_rows(levins_results) |>
  mutate(var = labels[var])
levins_order <- levins_df |>
  filter(ecotype == "LowT") |>
  arrange(B_std) |>
  pull(var)
levins_df$var <- factor(levins_df$var, levels = levins_order)

p_b <- ggplot(levins_df, aes(B_std, var, fill = ecotype)) +
  geom_col(position = position_dodge(0.68), width = 0.58) +
  scale_fill_manual(values = c(LowT = col_lowt, HighT = col_hight)) +
  labs(x = "Standardized Levins' B", y = NULL) +
  base_theme +
  theme(legend.position = "top", legend.key.width = unit(8, "pt"))

overlap_results <- lapply(env_vars, function(v) {
  breaks <- seq(min(dat[[v]]), max(dat[[v]]), length.out = 21)
  p_l <- tapply(lowt$AOA_ab, cut(lowt[[v]], breaks, include.lowest = TRUE), sum, default = 0)
  p_h <- tapply(hight$AOA_ab, cut(hight[[v]], breaks, include.lowest = TRUE), sum, default = 0)
  p_l[is.na(p_l)] <- 0
  p_h[is.na(p_h)] <- 0
  p_l <- p_l / sum(p_l)
  p_h <- p_h / sum(p_h)
  data.frame(var = labels[v], Schoener_D = 1 - 0.5 * sum(abs(p_l - p_h)))
})
overlap_df <- bind_rows(overlap_results) |>
  arrange(desc(Schoener_D))
overlap_df$var <- factor(overlap_df$var, levels = overlap_df$var)

p_c <- ggplot(overlap_df, aes(Schoener_D, var)) +
  geom_col(width = 0.58, fill = "#8E5BB7") +
  geom_text(aes(label = sprintf("%.3f", Schoener_D)),
            hjust = 1.1, size = 2.4, colour = "white") +
  scale_x_continuous(limits = c(0, 0.72), breaks = seq(0, 0.6, 0.2),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Schoener's D", y = NULL) +
  base_theme + theme(legend.position = "none")

supp_fig <- p_a + p_b + p_c +
  plot_layout(widths = c(1.25, 1.08, 0.95)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 10))
  )

ggsave(
  file.path(out_dir, "Supplementary_Fig_S3_thermal_niche.pdf"),
  supp_fig, width = 190, height = 65, units = "mm",
  device = cairo_pdf, bg = "white"
)

ggsave(
  file.path(out_dir, ".Supplementary_Fig_S3_thermal_niche_preview.png"),
  supp_fig, width = 190, height = 65, units = "mm", dpi = 200, bg = "white"
)
