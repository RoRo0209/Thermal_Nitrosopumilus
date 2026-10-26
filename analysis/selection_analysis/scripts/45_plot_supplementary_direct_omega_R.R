# V4 figure numbering is documented in docs/FIGURE_MAP.md; legacy filenames are retained.
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", winslash = "/", mustWork = TRUE)


suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

indir <- file.path(root, "selection_analysis/results/39_direct_hight_lowt_audit")
outdir <- file.path(root, "selection_analysis/results/44_supplementary_direct_omega_R")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

all <- read.delim(file.path(indir, "direct_omega_with_panelc_comparison.tsv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
hits <- read.delim(file.path(indir, "direct_omega_25_stable_FDR_hits.tsv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
excluded <- read.delim(
  file.path(indir, "direct_omega_3_FDR_hits_excluded_for_convergence.tsv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
stopifnot(nrow(all) == 110L, nrow(hits) == 25L, nrow(excluded) == 3L)

all$direct_log2_ratio <- log2(all$omega_HighT / all$omega_LowT)
all$minus_log10_fdr <- -log10(pmax(all$BH_FDR, .Machine$double.xmin))
all$status <- "Other genes"
all$status[all$split_family %in%
             hits$split_family[hits$direction == "HighT_lower"]] <-
  "Stable HighT-lower"
all$status[all$split_family == "group_2732"] <- "group_2732"
all$status[all$split_family %in% excluded$split_family] <-
  "FDR-positive, convergence-excluded"
stopifnot(sum(all$status == "Stable HighT-lower") == 24L)
stopifnot(sum(all$status == "group_2732") == 1L)
stopifnot(sum(all$status == "FDR-positive, convergence-excluded") == 3L)

status_levels <- c(
  "Other genes", "Stable HighT-lower", "group_2732",
  "FDR-positive, convergence-excluded"
)
all$status <- factor(all$status, levels = status_levels)

label_map <- c(
  eif6 = "eif6", pgk = "pgk", sufS_1 = "sufS", csl4 = "csl4",
  rpoD = "rpoD", rps9 = "rps9", group_2732 = "group_2732"
)
lab <- all[all$split_family %in% names(label_map), ]
lab$label <- unname(label_map[lab$split_family])
stopifnot(nrow(lab) == 7L)

cols <- c(
  "Other genes" = "#D1D1D1",
  "Stable HighT-lower" = "#66508C",
  "group_2732" = "#A34A5A",
  "FDR-positive, convergence-excluded" = "#30283D"
)
fills <- c(cols[1:3], "FDR-positive, convergence-excluded" = "white")
shapes <- c(
  "Other genes" = 16, "Stable HighT-lower" = 16,
  "group_2732" = 18, "FDR-positive, convergence-excluded" = 21
)
sizes <- c(
  "Other genes" = 1.55, "Stable HighT-lower" = 2.7,
  "group_2732" = 3.0, "FDR-positive, convergence-excluded" = 2.8
)

theme_pub <- theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.35, colour = "#333333"),
    axis.text = element_text(size = 6.7, colour = "#30283D"),
    axis.title = element_text(size = 7.2, colour = "#30283D"),
    legend.title = element_blank(),
    legend.text = element_text(size = 5.8),
    legend.key.width = grid::unit(8, "pt"),
    legend.spacing.x = grid::unit(2, "pt"),
    plot.tag = element_text(size = 8.5, face = "bold"),
    plot.margin = margin(5, 6, 3, 5)
  )

pa <- ggplot(
  all,
  aes(direct_log2_ratio, minus_log10_fdr,
      colour = status, fill = status, shape = status, size = status)
) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = 0.35, colour = "#777777") +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#555555") +
  geom_point(stroke = 0.75) +
  geom_text(
    data = lab, aes(label = label), inherit.aes = FALSE,
    x = lab$direct_log2_ratio, y = lab$minus_log10_fdr + 0.32,
    size = 2.35, fontface = ifelse(lab$split_family == "group_2732", "plain", "italic"),
    colour = "#30283D", check_overlap = FALSE
  ) +
  scale_colour_manual(values = cols, drop = FALSE) +
  scale_fill_manual(values = fills, drop = FALSE) +
  scale_shape_manual(values = shapes, drop = FALSE) +
  scale_size_manual(values = sizes, drop = FALSE) +
  guides(
    colour = guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(size = c(1.8, 2.5, 2.7, 2.5))
    ),
    fill = "none", shape = "none", size = "none"
  ) +
  labs(
    x = expression(log[2](omega[HighT] / omega[LowT])),
    y = expression(-log[10]("BH-FDR"))
  ) +
  theme_pub

pearson_r <- cor(all$panelc_log2_ratio, all$direct_log2_ratio, method = "pearson")
spearman_rho <- cor(all$panelc_log2_ratio, all$direct_log2_ratio, method = "spearman")
direction_match <- sum(sign(all$panelc_log2_ratio) == sign(all$direct_log2_ratio))
stopifnot(direction_match == 110L)
stat_label <- sprintf(
  "Pearson r = %.3f\nSpearman rho = %.3f\n110/110 direction consistent",
  pearson_r, spearman_rho
)

pb <- ggplot(all, aes(panelc_log2_ratio, direct_log2_ratio)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.45, colour = "#777777") +
  geom_point(size = 1.8, colour = "#66508C", alpha = 0.78) +
  annotate(
    "label", x = -Inf, y = Inf, label = stat_label,
    hjust = -0.08, vjust = 1.12, size = 2.35,
    family = "Arial", colour = "#30283D",
    fill = "white", linewidth = 0
  ) +
  coord_equal() +
  labs(
    x = expression("Old-model " * log[2](omega[HighT] / omega[LowT])),
    y = expression("Direct-model " * log[2](omega[HighT] / omega[LowT]))
  ) +
  theme_pub +
  theme(legend.position = "none")

fig <- (pa | pb) +
  plot_layout(widths = c(1.12, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

stem <- file.path(outdir, "Supplementary_Figure_direct_omega")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 183 / 25.4,
                     height = 92 / 25.4, family = "Arial", pointsize = 8)
print(fig)
grDevices::dev.off()
grDevices::png(paste0(stem, ".png"), width = 183 / 25.4,
               height = 92 / 25.4, units = "in", res = 600,
               type = "cairo", bg = "white")
print(fig)
grDevices::dev.off()
