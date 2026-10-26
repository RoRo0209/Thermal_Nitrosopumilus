#!/usr/bin/env Rscript

# Figure 7c: ranked HighT/LowT foreground omega ratios with six concise,
# functionally supported direct labels. Drawing and export are R-only.

project_root <- normalizePath(
  if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[1] else ".",
  winslash = "/", mustWork = TRUE
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

outdir <- file.path(
  project_root,
  "selection_analysis/results/34_panelc_annotated_R"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

source_file <- file.path(
  project_root,
  "selection_analysis/results/16_dnds_lineage_history/source_data_branch_model.tsv"
)
annotation_file <- file.path(
  project_root,
  "selection_analysis/results/32_blastkoala_integration/",
  "Figure7c_110_families_Kofam_BlastKOALA.tsv"
)
direct_hit_file <- file.path(
  project_root,
  "selection_analysis/results/39_direct_hight_lowt_audit/",
  "direct_omega_25_stable_FDR_hits.tsv"
)
direct_excluded_file <- file.path(
  project_root,
  "selection_analysis/results/39_direct_hight_lowt_audit/",
  "direct_omega_3_FDR_hits_excluded_for_convergence.tsv"
)

COL_POINT <- "#66508C"
COL_DARK <- "#30283D"
COL_MUTED <- "#696969"
COL_GRID <- "#E8E8E8"
COL_BELOW <- "#F5F1FA"
COL_ABOVE <- "#F5F5F5"
COL_SEGMENT <- "#9586AD"

raw <- read.delim(source_file, check.names = FALSE, stringsAsFactors = FALSE)
ann <- read.delim(
  annotation_file, check.names = FALSE, stringsAsFactors = FALSE,
  quote = "", comment.char = ""
)
direct_hits <- read.delim(
  direct_hit_file, check.names = FALSE, stringsAsFactors = FALSE
)
direct_excluded <- read.delim(
  direct_excluded_file, check.names = FALSE, stringsAsFactors = FALSE
)
if (nrow(direct_hits) != 25L || nrow(direct_excluded) != 3L) {
  stop("Expected 25 stable direct hits and 3 convergence-excluded hits.")
}

needed <- c("split_family", "annotation", "ecotype", "omega_foreground")
missing_cols <- setdiff(needed, names(raw))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

high <- raw[
  raw$ecotype == "HighT",
  c("split_family", "annotation", "omega_foreground")
]
low <- raw[
  raw$ecotype == "LowT",
  c("split_family", "omega_foreground")
]
names(high)[3] <- "omega_HighT"
names(low)[2] <- "omega_LowT"

wide <- merge(high, low, by = "split_family", all = FALSE, sort = FALSE)
if (nrow(wide) != 110L) {
  stop("Expected 110 tested families; observed ", nrow(wide))
}

wide$omega_ratio_HighT_LowT <- wide$omega_HighT / wide$omega_LowT
wide$effect_log2 <- log2(wide$omega_ratio_HighT_LowT)
wide <- wide[order(wide$effect_log2, wide$split_family), ]
wide$rank <- seq_len(nrow(wide))

annotation_columns <- c(
  "split_family", "Integrated_KO", "Recommended_annotation",
  "panelC_direct_support", "HighT_branch_FDR", "LowT_branch_FDR"
)
if (!all(annotation_columns %in% names(ann))) {
  stop("Integrated annotation table lacks required columns.")
}
wide <- merge(
  wide,
  ann[, annotation_columns],
  by = "split_family",
  all.x = TRUE,
  sort = FALSE
)
wide <- wide[order(wide$rank), ]
wide$direct_status <- "Other"
wide$direct_status[wide$split_family %in% direct_hits$split_family] <-
  "Stable accepted FDR hit"
wide$direct_status[wide$split_family %in% direct_excluded$split_family] <-
  "FDR-positive, convergence-excluded"
if (sum(wide$direct_status == "Stable accepted FDR hit") != 25L ||
    sum(wide$direct_status == "FDR-positive, convergence-excluded") != 3L) {
  stop("Direct-model classifications do not map uniquely onto the 110 genes.")
}

n_h_lower <- sum(wide$omega_ratio_HighT_LowT < 1)
n_l_lower <- sum(wide$omega_ratio_HighT_LowT > 1)
if (!identical(c(n_h_lower, n_l_lower), c(81L, 29L))) {
  stop("Unexpected direction counts: ", n_h_lower, " and ", n_l_lower)
}

median_ratio <- median(wide$omega_ratio_HighT_LowT)
median_log2 <- log2(median_ratio)
sign_p <- binom.test(
  n_h_lower, nrow(wide), p = 0.5, alternative = "two.sided"
)$p.value

label_spec <- data.frame(
  split_family = c("eif6", "sufS_1", "csl4", "pgk", "rps9", "rpoD"),
  display_label = c("eif6", "sufS", "csl4", "pgk", "rps9", "rpoD"),
  label_x = c(2.7, 12.5, 4.4, 10.2, 10.0, 29.0),
  label_y = c(-2.67, -2.18, -1.31, -0.88, -2.87, -0.86),
  hjust = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  stringsAsFactors = FALSE
)

label_df <- merge(
  wide,
  label_spec,
  by = "split_family",
  all = FALSE,
  sort = FALSE
)
label_df <- label_df[match(label_spec$split_family, label_df$split_family), ]

if (nrow(label_df) != 6L ||
    any(!label_df$split_family %in% direct_hits$split_family)) {
  stop("The six labels must all be stable accepted direct-model hits.")
}

wide$display_label <- ""
wide$display_label[match(label_df$split_family, wide$split_family)] <-
  label_df$display_label

ratio_breaks <- c(0.125, 0.25, 0.5, 1, 2)

p <- ggplot(wide, aes(x = rank, y = effect_log2)) +
  annotate(
    "rect", xmin = 0.5, xmax = 110.5, ymin = -Inf, ymax = 0,
    fill = COL_BELOW
  ) +
  annotate(
    "rect", xmin = 0.5, xmax = 110.5, ymin = 0, ymax = Inf,
    fill = COL_ABOVE
  ) +
  geom_hline(yintercept = 0, colour = "#555555", linewidth = 0.55) +
  geom_segment(
    aes(xend = rank, y = 0, yend = effect_log2),
    colour = alpha(COL_POINT, 0.46), linewidth = 0.32
  ) +
  geom_point(
    shape = 21, size = 1.62, stroke = 0.24,
    fill = alpha(COL_POINT, 0.82), colour = "white"
  ) +
  geom_point(
    data = wide[wide$direct_status == "Stable accepted FDR hit", ],
    aes(x = rank, y = effect_log2),
    inherit.aes = FALSE,
    shape = 21, size = 2.22, stroke = 0.34,
    fill = COL_POINT, colour = "white"
  ) +
  geom_point(
    data = wide[
      wide$direct_status == "FDR-positive, convergence-excluded", ],
    aes(x = rank, y = effect_log2),
    inherit.aes = FALSE,
    shape = 21, size = 2.35, stroke = 0.68,
    fill = "white", colour = COL_DARK
  ) +
  geom_hline(
    yintercept = median_log2, linetype = "dashed",
    colour = COL_DARK, linewidth = 0.46
  ) +
  annotate(
    "label", x = 47, y = median_log2 + 0.06,
    label = "median",
    hjust = 0, vjust = 0, size = 2.25, family = "Arial",
    colour = COL_DARK, fill = alpha("white", 0.90), linewidth = 0,
    label.padding = grid::unit(0.05, "lines")
  ) +
  geom_segment(
    data = label_df,
    aes(x = rank, y = effect_log2, xend = label_x, yend = label_y),
    inherit.aes = FALSE,
    colour = COL_SEGMENT, linewidth = 0.32, lineend = "round"
  ) +
  geom_point(
    data = label_df,
    aes(x = rank, y = effect_log2),
    inherit.aes = FALSE,
    shape = 21, size = 2.15, stroke = 0.35,
    fill = COL_DARK, colour = "white"
  ) +
  geom_text(
    data = label_df,
    aes(x = label_x, y = label_y, label = display_label, hjust = hjust),
    inherit.aes = FALSE,
    family = "Arial", fontface = "italic",
    size = 2.45, colour = COL_DARK, vjust = 0.5
  ) +
  scale_x_continuous(
    limits = c(0.5, 110.5),
    breaks = c(1, 25, 50, 75, 100, 110),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-3.20, 1.35),
    breaks = log2(ratio_breaks),
    labels = c("0.125", "0.25", "0.5", "1", "2"),
    expand = c(0, 0)
  ) +
  labs(
    title = expression("c   Core-gene foreground " * omega * " ratios"),
    subtitle = paste0(
      "81/110 < 1; median = 0.69; sign test  P = 7.3e-7\n",
      "Direct H\u2013L test: 25 stable FDR hits; 24 HighT-lower"
    ),
    x = expression("Genes ranked by HighT / LowT " * omega * " ratio"),
    y = expression("HighT / LowT foreground " * omega * " ratio")
  ) +
  theme_minimal(base_size = 7.5, base_family = "Arial") +
  theme(
    text = element_text(colour = COL_DARK, family = "Arial"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = COL_GRID, linewidth = 0.28),
    axis.line = element_line(colour = "#555555", linewidth = 0.35),
    axis.ticks = element_line(colour = "#555555", linewidth = 0.35),
    axis.title = element_text(size = 7.5),
    axis.text = element_text(size = 7.0, colour = COL_DARK),
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 6.8, colour = COL_MUTED),
    plot.margin = margin(3, 5, 3, 4)
  )

save_pub <- function(plot, stem, width_mm = 183, height_mm = 62) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = "Arial", pointsize = 8, onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()

  grDevices::png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300, type = "cairo", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
}

stem <- file.path(
  outdir, "Figure7c_ranked_omega_ratio_annotated_R_directValidated"
)
save_pub(p, stem)
cat("Wrote direct-validated Figure 7c PDF and PNG to:", outdir, "\n")
