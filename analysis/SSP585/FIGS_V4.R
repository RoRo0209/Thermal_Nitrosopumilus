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

## =========================================================
## FIGS_V4.R — Temperature–Chlorophyll environmental space
## 2D density plot with ALE-informed quadrant partitions
## =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

## ---------------------------------------------------------
## 0) Paths
## ---------------------------------------------------------
base <- thermal_path("SSP585")

high_map <- file.path(base, "NP_highT", "AOA_future_projection_outputs",
                      "GLOBAL_NEARSHORE_MAPS")
tag <- "nearshore200km_KEEPCHL_SSP585_2080-2100"

## ---------------------------------------------------------
## 1) Load data
## ---------------------------------------------------------
h_delta <- read_tsv(file.path(high_map, paste0("GLOBAL_nearshore_pred_delta_", tag, ".tsv")),
                    show_col_types = FALSE)
h_cur   <- read_tsv(file.path(high_map, paste0("GLOBAL_nearshore_pred_current_", tag, ".tsv")),
                    show_col_types = FALSE)
h_fut   <- read_tsv(file.path(high_map, paste0("GLOBAL_nearshore_pred_future_",  tag, ".tsv")),
                    show_col_types = FALSE)

# Fail rather than combine event classes and environmental values from different cells.
thermal_coords <- c("Longitude", "Latitude")
for (x in list(h_fut)) {
  if (!all(thermal_coords %in% names(x)) || !all(thermal_coords %in% names(h_delta)))
    stop("All prediction tables require Longitude and Latitude columns.")
  if (anyDuplicated(x[thermal_coords]) || anyDuplicated(h_delta[thermal_coords]) ||
      !identical(as.data.frame(x[thermal_coords]), as.data.frame(h_delta[thermal_coords])))
    stop("Prediction rows do not match by unique coordinates; regenerate aligned tables before plotting.")
}

df <- tibble(
  delta_highT = as.numeric(h_delta$delta_future_minus_current),
  sst_fut = as.numeric(h_fut$temp),
  chl_fut = as.numeric(h_fut$chl)
) %>%
  filter(is.finite(delta_highT), is.finite(sst_fut), is.finite(chl_fut)) %>%
  mutate(
    change_class = factor(
      case_when(delta_highT >  0.2 ~ "Increase",
                delta_highT < -0.2 ~ "Decrease",
                TRUE               ~ "Stable"),
      levels = c("Increase", "Stable", "Decrease"))
  )

message("Grid cells: ", nrow(df))

## ---------------------------------------------------------
## 2) ALE-informed boundaries
##    Temperature: 22°C (HighT ALE steep rise begins ~22°C)
##    Chlorophyll: 1.5 mg/m³ (HighT ALE steep rise above ~1.5)
## ---------------------------------------------------------
sst_cut <- 22
chl_cut <- 1.5

## Quadrant labels
quad_label <- function(sst, chl) {
  case_when(
    sst >= sst_cut & chl >= chl_cut ~ "Warm + Productive",
    sst >= sst_cut & chl <  chl_cut ~ "Warm + Low-productivity",
    sst <  sst_cut & chl >= chl_cut ~ "Cool + Productive",
    sst <  sst_cut & chl <  chl_cut ~ "Cool + Low-productivity"
  )
}

df <- df %>%
  mutate(quadrant = factor(quad_label(sst_fut, chl_fut),
                           levels = c("Warm + Productive",
                                      "Warm + Low-productivity",
                                      "Cool + Productive",
                                      "Cool + Low-productivity")))

## ---------------------------------------------------------
## 3) Quadrant statistics table
## ---------------------------------------------------------
quad_tbl <- df %>%
  group_by(change_class, quadrant) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(change_class) %>%
  mutate(fraction = n / sum(n)) %>%
  ungroup() %>%
  arrange(change_class, quadrant)

message("\n=== Quadrant occupancy ===")
print(as.data.frame(quad_tbl), row.names = FALSE)

write_tsv(quad_tbl, file.path(base, "FIGS_V4_quadrant_table.tsv"))

## ---------------------------------------------------------
## 4) Panel A: Faceted 2D density fill by response group
## ---------------------------------------------------------
change_colors <- c(Increase = "#D73027", Stable = "grey50", Decrease = "#4575B4")

## Clip chl for better visualization
chl_upper <- quantile(df$chl_fut, 0.98)

## Shared axis limits
xlims <- c(-2, 38)
ylims <- c(0, chl_upper)

## Facet labels with n
facet_labels <- df %>%
  group_by(change_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(change_class, "  (n = ", format(n, big.mark = ","), ")"))
facet_named <- setNames(facet_labels$label, facet_labels$change_class)

## Quadrant label positions
quad_annot <- tibble(
  x = c(36, 36, 2, 2),
  y = c(chl_upper * 0.92, ylims[1] + 0.2, chl_upper * 0.92, ylims[1] + 0.2),
  label = c("Warm +\nProductive", "Warm +\nLow-productivity",
            "Cool +\nProductive", "Cool +\nLow-productivity")
)

p_main <- ggplot(df, aes(x = sst_fut, y = chl_fut)) +
  ## 2D filled density (no legend — density not comparable across facets)
  stat_density_2d_filled(contour_var = "ndensity", bins = 8, alpha = 0.8) +
  ## ALE-informed boundary lines
  geom_vline(xintercept = sst_cut, linetype = "dashed",
             color = "grey20", linewidth = 0.5) +
  geom_hline(yintercept = chl_cut, linetype = "dashed",
             color = "grey20", linewidth = 0.5) +
  ## Quadrant labels
  geom_text(data = quad_annot, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 2.3, color = "grey30",
            fontface = "italic", lineheight = 0.85) +
  facet_wrap(~ change_class, nrow = 1,
             labeller = labeller(change_class = facet_named)) +
  scale_fill_manual(
    values = colorRampPalette(c("white", "#FEE08B", "#FC8D59",
                                "#D73027", "#67001F"))(9),
    guide = "none") +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  labs(x = expression("Future SST (\u00B0C)"),
       y = expression("Future Chl (mg m"^{-3}*")")) +
  theme_classic(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 10.5),
        strip.background = element_rect(fill = "grey95", color = NA),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        panel.grid.major = element_line(color = "grey95", linewidth = 0.2),
        panel.spacing = unit(0.8, "lines"))

## ---------------------------------------------------------
## 5) Panel B: Quadrant fraction stacked bar chart
## ---------------------------------------------------------
## Reorder quadrant factor for consistent display
quad_order <- c("Cool + Low-productivity",
                "Cool + Productive",
                "Warm + Low-productivity",
                "Warm + Productive")

quad_tbl <- quad_tbl %>%
  mutate(quadrant = factor(quadrant, levels = quad_order))

quad_colors <- c("Cool + Low-productivity"  = "#ABD9E9",
                 "Cool + Productive"        = "#74ADD1",
                 "Warm + Low-productivity"  = "#FDAE61",
                 "Warm + Productive"        = "#D73027")

## Add percentage labels inside bars
quad_tbl <- quad_tbl %>%
  group_by(change_class) %>%
  arrange(change_class, quadrant) %>%
  mutate(cum_frac = cumsum(fraction),
         mid_frac = cum_frac - fraction / 2,
         pct_label = ifelse(fraction >= 0.05,
                            paste0(round(fraction * 100, 0), "%"), "")) %>%
  ungroup()

p_bar <- ggplot(quad_tbl, aes(x = change_class, y = fraction, fill = quadrant)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(y = mid_frac, label = pct_label),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = quad_colors, name = "Environmental sector",
                    breaks = rev(quad_order)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "HighT abundance change",
       y = "Fraction of grid cells") +
  theme_classic(base_size = 10) +
  theme(legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9.5, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

## ---------------------------------------------------------
## 6) Combine and save
## ---------------------------------------------------------
combined <- p_main / p_bar +
  plot_layout(heights = c(2, 1.2)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(", tag_suffix = ")",
                  theme = theme(plot.tag = element_text(face = "bold", size = 13)))

out_pdf <- file.path(base, "FIGS_V4_env_space_v2.pdf")

ggsave(out_pdf, combined,
       width = 13, height = 10, dpi = 300,
       device = cairo_pdf, bg = "white")

message("\nFigure saved: ", out_pdf)

## ---------------------------------------------------------
## 7) Results text and figure legend
## ---------------------------------------------------------
q_inc_wp <- quad_tbl %>% filter(change_class == "Increase", quadrant == "Warm + Productive")
q_inc_wl <- quad_tbl %>% filter(change_class == "Increase", quadrant == "Warm + Low-productivity")
q_dec_wp <- quad_tbl %>% filter(change_class == "Decrease", quadrant == "Warm + Productive")
q_dec_wl <- quad_tbl %>% filter(change_class == "Decrease", quadrant == "Warm + Low-productivity")
q_inc_cp <- quad_tbl %>% filter(change_class == "Increase", quadrant == "Cool + Productive")

legend_txt <- paste0(
  "\n=== Figure legend ===\n\n",
  "Fig. SXX. Future environmental position by HighT response group. ",
  "(a) Faceted density plots showing the distribution of nearshore grid cells in ",
  "future SST\u2013chlorophyll space, stratified by projected HighT abundance change ",
  "(Increase: \u0394log\u2081\u2080 > 0.2; Stable: |\u0394log\u2081\u2080| \u2264 0.2; ",
  "Decrease: \u0394log\u2081\u2080 < \u22120.2). Dashed lines indicate ALE-informed ",
  "response-partition boundaries (SST = ", sst_cut, "\u00B0C; Chl = ", chl_cut,
  " mg m\u207B\u00B3) that divide the environmental space into four sectors. ",
  "(b) Stacked bar chart summarising the fraction of grid cells in each ",
  "environmental sector by response group.\n"
)

results_txt <- paste0(
  "\n=== Results text ===\n\n",
  "HighT increases were more often associated with productive sectors of ",
  "environmental space: ", round(q_inc_wp$fraction * 100, 1),
  "% of Increase cells fell within the warm, productive sector and ",
  round(q_inc_cp$fraction * 100, 1),
  "% within the cool, productive sector, whereas ",
  round(q_dec_wl$fraction * 100, 1),
  "% of Decrease cells were concentrated in warm but low-productivity environments. ",
  "These results indicate that the future response of the HighT ecotype depended on ",
  "final position in joint SST\u2013chlorophyll space rather than on temperature change alone, ",
  "explaining why warming was insufficient to drive net HighT expansion.\n"
)

cat(legend_txt)
cat(results_txt)
writeLines(c(legend_txt, results_txt), file.path(base, "FIGS_V4_results_text.txt"))

message("\nDONE")
