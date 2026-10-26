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
setwd(thermal_path("genome_character/PCA-AA"))

## ============================================================
## Amino-acid usage PCA: points + convex hull + biplot arrows
## - Arrows colored by Contribution (%) with colorbar (Fig B-like)
## Output: AA_PCA_biplot_contrib.pdf / .png
## ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

## =========================
## 1) Input files (按需改路径)
## =========================
aa_file   <- thermal_required_env("THERMAL_AA_USAGE")
meta_file <- thermal_required_env("THERMAL_AA_METADATA")

## =========================
## 2) Read data
## =========================
aa <- read_tsv(aa_file, show_col_types = FALSE)
meta <- read_excel(meta_file)

stopifnot(ncol(aa) >= 3)
stopifnot(all(c("FileName","group") %in% colnames(meta)))

dat <- aa %>%
  left_join(meta %>% transmute(`Genome ID` = FileName, group), by = "Genome ID") %>%
  mutate(group = ifelse(is.na(group), "Other", as.character(group)))

## =========================
## 3) PCA
## =========================
aa_cols <- setdiff(colnames(dat), c("Genome ID", "group"))
dat[aa_cols] <- lapply(dat[aa_cols], function(x) readr::parse_number(as.character(x)))

pca <- prcomp(dat[, aa_cols], center = TRUE, scale. = TRUE)

var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
pc1_lab <- sprintf("PC1 (%.1f%%)", 100 * var_exp[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", 100 * var_exp[2])

scores <- as.data.frame(pca$x[, 1:2]) %>%
  mutate(`Genome ID` = dat$`Genome ID`,
         group = dat$group)

## =========================
## 4) Convex hull per group
## =========================
hulls <- scores %>%
  group_by(group) %>%
  filter(n() >= 3) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

## =========================
## 5) Loadings + contributions (for coloring arrows)
## =========================
loadings <- as.data.frame(pca$rotation[, 1:2])
loadings$AA <- rownames(loadings)
colnames(loadings)[1:2] <- c("PC1", "PC2")

contrib <- loadings %>%
  mutate(
    PC1_contrib = (PC1^2) / sum(PC1^2) * 100,
    PC2_contrib = (PC2^2) / sum(PC2^2) * 100,
    ## 一根箭头一个颜色：取 PC1/PC2 中更大的贡献
    contrib_max = pmax(PC1_contrib, PC2_contrib)
  )

## 选“总体贡献”最大的若干 AA 来画箭头（避免太挤）
top_n <- 12   # 你可以改成 20 画全套AA
load_plot <- contrib %>%
  arrange(desc(contrib_max)) %>%
  slice(1:top_n)

message("Top AA (by max contribution on PC1/PC2):")
print(load_plot %>% select(AA, PC1_contrib, PC2_contrib, contrib_max) %>%
        mutate(across(where(is.numeric), ~round(.x, 1))))

## =========================
## 6) Arrow scaling (auto)
## =========================
range_scores <- max(abs(c(scores$PC1, scores$PC2)), na.rm = TRUE)
range_load   <- max(sqrt(loadings$PC1^2 + loadings$PC2^2), na.rm = TRUE)
arrow_scale  <- (range_scores * 0.35) / range_load

## =========================
## 7) Group colors
## =========================
custom_cols <- c(
  "cold-wide" = "#8EC6E8",
  "warm-wide" = "#F4A261",
  "WCA"       = "#9467bd",
  "Other"     = "grey60"
)

## =========================
## 8) Plot: hull + points + arrows colored by contribution
## =========================
p <- ggplot() +
  geom_polygon(
    data = hulls,
    aes(x = PC1, y = PC2, fill = group, color = group, group = group),
    alpha = 0.20,
    linewidth = 0.7
  ) +
  geom_point(
    data = scores,
    aes(x = PC1, y = PC2, color = group),
    size = 2.1,
    alpha = 0.9
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  
  ## arrows: map color to contrib_max (with colorbar)
  geom_segment(
    data = load_plot,
    aes(
      x = 0, y = 0,
      xend = PC1 * arrow_scale,
      yend = PC2 * arrow_scale,
      color = contrib_max
    ),
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = load_plot,
    aes(
      x = PC1 * arrow_scale * 1.10,
      y = PC2 * arrow_scale * 1.10,
      label = AA
    ),
    size = 4,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  ) +
  
  ## group colors (for points/hulls) — use new scales trick:
  ## easiest: set group colors first, then override arrow color scale via ggnewscale
  ## BUT we avoid extra package by doing this:
  ## - draw points/hulls with fixed colors using scale_*_manual
  ## - arrows use a separate aesthetic mapping "color = contrib_max" already,
  ##   so group scale_color_manual cannot be used directly.
  ## Solution: set point colors via "colour" outside mapping:
  NULL

## ---- rebuild p with safe two-color-scale approach (no ggnewscale) ----
## We will draw points/hulls with "group_col" computed column (fixed),
## then use color scale only for arrows.

scores2 <- scores %>% mutate(group_col = custom_cols[as.character(group)])
hulls2  <- hulls  %>% mutate(group_col = custom_cols[as.character(group)])

p <- ggplot() +
  geom_polygon(
    data = hulls2,
    aes(x = PC1, y = PC2, group = group),
    fill = hulls2$group_col,
    color = hulls2$group_col,
    alpha = 0.20,
    linewidth = 0.7
  ) +
  geom_point(
    data = scores2,
    aes(x = PC1, y = PC2),
    color = scores2$group_col,
    size = 2.1,
    alpha = 0.9
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  
  geom_segment(
    data = load_plot,
    aes(
      x = 0, y = 0,
      xend = PC1 * arrow_scale,
      yend = PC2 * arrow_scale,
      color = contrib_max
    ),
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.9
  ) +
  geom_text(
    data = load_plot,
    aes(
      x = PC1 * arrow_scale * 1.10,
      y = PC2 * arrow_scale * 1.10,
      label = AA
    ),
    size = 4,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  ) +
  
  ## contribution colorbar (Fig B-like)
  scale_color_gradient(
    name = "Contribution (%)",
    low  = "#2DC7D3",
    high = "#F05A28"
  ) +
  
  labs(x = pc1_lab, y = pc2_lab) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    legend.position = "right"
  )

print(p)

ggsave("AA_PCA_biplot_contrib.pdf", p, width = 6.6, height = 5.0)
ggsave("AA_PCA_biplot_contrib.png", p, width = 6.6, height = 5.0, dpi = 300)

write.csv(contrib %>% arrange(desc(contrib_max)),
          "AA_PCA_loading_contributions.csv", row.names = FALSE)

## Done.

## ==========
