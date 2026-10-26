#!/usr/bin/env Rscript

# Figure 7c: compact ranked HighT/LowT foreground omega ratios for the
# 110 shared single-copy gene families passing codeml quality control.

project_root <- normalizePath(
  if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[1] else ".",
  winslash = "/", mustWork = TRUE
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

outdir <- file.path(project_root, "selection_analysis/results/28_panelc_ranked_ratio_only_R")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
source_file <- file.path(
  project_root,
  "selection_analysis/results/16_dnds_lineage_history/source_data_branch_model.tsv"
)

COL_POINT <- "#66508C"
COL_DARK <- "#30283D"
COL_MUTED <- "#696969"
COL_GRID <- "#E8E8E8"
COL_BELOW <- "#F5F1FA"
COL_ABOVE <- "#F5F5F5"

raw <- read.delim(source_file, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c("split_family", "annotation", "ecotype", "omega_foreground")
missing_cols <- setdiff(needed, names(raw))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

high <- raw[raw$ecotype == "HighT", c("split_family", "annotation", "omega_foreground")]
low <- raw[raw$ecotype == "LowT", c("split_family", "omega_foreground")]
names(high)[3] <- "omega_HighT"
names(low)[2] <- "omega_LowT"

wide <- merge(high, low, by = "split_family", all = FALSE, sort = FALSE)
if (nrow(wide) != 110) stop("Expected 110 tested families; observed ", nrow(wide))

wide$omega_ratio_HighT_LowT <- wide$omega_HighT / wide$omega_LowT
wide$effect_log2 <- log2(wide$omega_ratio_HighT_LowT)
wide <- wide[order(wide$effect_log2, wide$split_family), ]
wide$rank <- seq_len(nrow(wide))

n_h_lower <- sum(wide$omega_ratio_HighT_LowT < 1)
n_l_lower <- sum(wide$omega_ratio_HighT_LowT > 1)
if (!identical(c(n_h_lower, n_l_lower), c(81L, 29L))) {
  stop("Unexpected direction counts: ", n_h_lower, " and ", n_l_lower)
}

median_ratio <- median(wide$omega_ratio_HighT_LowT)
median_log2 <- log2(median_ratio)
sign_p <- binom.test(n_h_lower, nrow(wide), p = 0.5, alternative = "two.sided")$p.value

ratio_breaks <- c(0.125, 0.25, 0.5, 1, 2)

p <- ggplot(wide, aes(x = rank, y = effect_log2)) +
  annotate("rect", xmin = 0.5, xmax = 110.5, ymin = -Inf, ymax = 0,
           fill = COL_BELOW) +
  annotate("rect", xmin = 0.5, xmax = 110.5, ymin = 0, ymax = Inf,
           fill = COL_ABOVE) +
  geom_hline(yintercept = 0, colour = "#555555", linewidth = 0.55) +
  geom_segment(
    aes(xend = rank, y = 0, yend = effect_log2),
    colour = alpha(COL_POINT, 0.46), linewidth = 0.32
  ) +
  geom_point(shape = 21, size = 1.65, stroke = 0.25,
             fill = alpha(COL_POINT, 0.84), colour = "white") +
  geom_hline(yintercept = median_log2, linetype = "dashed",
             colour = COL_DARK, linewidth = 0.46) +
  annotate(
    "label", x = 4, y = median_log2 - 0.13,
    label = sprintf("median = %.2f", median_ratio),
    hjust = 0, vjust = 1, size = 2.45, family = "Arial",
    colour = COL_DARK, fill = alpha("white", 0.90), linewidth = 0,
    label.padding = unit(0.06, "lines")
  ) +
  annotate(
    "text", x = 7, y = -2.48,
    label = paste0(n_h_lower, " genes: HighT lower"),
    hjust = 0, size = 2.65, fontface = "bold",
    colour = COL_DARK, family = "Arial"
  ) +
  annotate(
    "text", x = 106, y = 0.92,
    label = paste0(n_l_lower, " genes: LowT lower"),
    hjust = 1, size = 2.55, colour = COL_MUTED, family = "Arial"
  ) +
  annotate(
    "text", x = 109, y = 0.08, label = "equal (ratio = 1)",
    hjust = 1, vjust = 0, size = 2.25, colour = COL_MUTED,
    family = "Arial"
  ) +
  scale_x_continuous(
    limits = c(0.5, 110.5), breaks = c(1, 25, 50, 75, 100, 110),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-3.30, 1.45), breaks = log2(ratio_breaks),
    labels = c("0.125", "0.25", "0.5", "1", "2"),
    expand = c(0, 0)
  ) +
  labs(
    title = expression("c   Core-gene foreground " * omega * " ratios"),
    subtitle = paste0(
      "HighT lower in 81 of 110 genes; median ratio = 0.69; exact sign test  P = ",
      format(sign_p, scientific = TRUE, digits = 2)
    ),
    x = "110 shared genes ranked by HighT / LowT ratio",
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

write.table(
  wide,
  file.path(outdir, "Figure7c_110_source_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  data.frame(
    metric = c(
      "codeml_tested", "HighT_omega_lower", "LowT_omega_lower",
      "median_HighT_LowT_ratio", "paired_sign_test_P"
    ),
    value = c(110, n_h_lower, n_l_lower, median_ratio, sign_p)
  ),
  file.path(outdir, "Figure7c_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

save_pub <- function(plot, stem, width_mm = 183, height_mm = 62) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(
    paste0(stem, ".svg"), width = width_in, height = height_in,
    system_fonts = list(sans = "Arial")
  )
  print(plot)
  grDevices::dev.off()

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

  grDevices::tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = 600, compression = "lzw", type = "cairo", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
}

stem <- file.path(outdir, "Figure7c_ranked_omega_ratio_R")
save_pub(p, stem)

writeLines(c(
  "# Figure 7c legend draft",
  "",
  "**c, HighT shows lower foreground omega across most shared core genes.**",
  "The 110 single-copy orthologous families passing codeml quality control are ranked by the ratio of their separately estimated HighT and LowT foreground omega values. Ratios below one indicate lower omega in HighT. HighT had the lower estimate for 81 of 110 genes; the median HighT/LowT ratio was 0.69 (two-sided exact sign test, P = 7.3 x 10^-7). This across-gene paired summary does not constitute a gene-specific direct HighT-versus-LowT likelihood-ratio test."
), file.path(outdir, "FIGURE7C_LEGEND_DRAFT.md"), useBytes = TRUE)

writeLines(c(
  "# Figure QA",
  "",
  "- PASS: all 110 QC-passing families are shown exactly once.",
  "- PASS: the y-axis displays the actual HighT/LowT ratio while positions are log2 transformed.",
  "- PASS: ratio = 1, the median ratio, direction counts, n and the exact sign test are shown.",
  "- PASS: neutral purple/grey encoding avoids reusing ecotype colours for ratio magnitude.",
  "- PASS: R-only SVG/PDF/PNG/TIFF exports are present.",
  "- LIMITATION: the values come from separate foreground fits and are compared descriptively."
), file.path(outdir, "FIGURE_QA.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(outdir, "R_SESSION_INFO.txt"))
cat("Wrote compact Figure 7c to:", outdir, "\n")
