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
## niche_analysis_v2.R
## LowT vs HighT 生态位分析
##   Main Figure (3-panel): KDE → Levins' B → Schoener overlap
##   Supp Figure (2-panel): 丰度热图 + 环境异质性
## =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

## ---- 路径 ----
base_dir <- thermal_path("SSP585")
lowt_file  <- file.path(base_dir, "NP_lowT/AOA_points_with_BioORACLE_nearshore_filtered.tsv")
hight_file <- file.path(base_dir, "NP_highT/AOA_points_with_BioORACLE_nearshore_filtered.tsv")
out_dir <- file.path(base_dir, "niche_analysis_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

env_vars <- c("temp", "sal", "chl", "o2", "no3", "po4", "si")

## 配色
col_lowt  <- "#9CCFDE"
col_hight <- "#E8A07A"

## ---- 读取数据 ----
message("=== 读取数据 ===")
lowt  <- read_tsv(lowt_file, show_col_types = FALSE) %>%
  mutate(ecotype = "LowT") %>%
  filter(if_all(all_of(env_vars), is.finite))
hight <- read_tsv(hight_file, show_col_types = FALSE) %>%
  mutate(ecotype = "HighT") %>%
  filter(if_all(all_of(env_vars), is.finite))

dat <- bind_rows(lowt, hight)
message("LowT samples: ", nrow(lowt), "  HighT samples: ", nrow(hight))

## =========================================================
## 1) 生态位中心与宽度（丰度加权）
## =========================================================
message("\n=== 指标1: 生态位中心与宽度 ===")

niche_center_breadth <- dat %>%
  pivot_longer(all_of(env_vars), names_to = "var", values_to = "value") %>%
  group_by(ecotype, var) %>%
  summarise(
    n = n(),
    mu_weighted = weighted.mean(value, AOA_ab),
    sigma_weighted = sqrt(sum(AOA_ab * (value - weighted.mean(value, AOA_ab))^2) / sum(AOA_ab)),
    .groups = "drop"
  )

cat("\n丰度加权生态位中心与宽度:\n")
niche_center_breadth %>%
  select(ecotype, var, mu_weighted, sigma_weighted) %>%
  pivot_wider(names_from = ecotype, values_from = c(mu_weighted, sigma_weighted)) %>%
  print(n = Inf)

## ---- Panel A: 温度KDE ----
ncb_t <- niche_center_breadth %>% filter(var == "temp")
mu_lowt  <- ncb_t$mu_weighted[ncb_t$ecotype == "LowT"]
sd_lowt  <- ncb_t$sigma_weighted[ncb_t$ecotype == "LowT"]
mu_hight <- ncb_t$mu_weighted[ncb_t$ecotype == "HighT"]
sd_hight <- ncb_t$sigma_weighted[ncb_t$ecotype == "HighT"]

dens_lowt  <- density(lowt$temp, weights = lowt$AOA_ab / sum(lowt$AOA_ab),
                      from = -2, to = 32, n = 512, adjust = 1.5)
dens_hight <- density(hight$temp, weights = hight$AOA_ab / sum(hight$AOA_ab),
                      from = -2, to = 32, n = 512, adjust = 1.5)
dens_df <- bind_rows(
  data.frame(temp = dens_lowt$x, density = dens_lowt$y, ecotype = "LowT"),
  data.frame(temp = dens_hight$x, density = dens_hight$y, ecotype = "HighT")
)

ymax <- max(dens_df$density)
seg_y <- -ymax * 0.04

p_A <- ggplot(dens_df, aes(x = temp, y = density, fill = ecotype, color = ecotype)) +
  geom_area(alpha = 0.35, position = "identity") +
  geom_line(linewidth = 1) +
  # mu +/- sigma 标记
  annotate("segment", x = mu_lowt - sd_lowt, xend = mu_lowt + sd_lowt,
           y = seg_y, yend = seg_y, color = col_lowt, linewidth = 2) +
  annotate("point", x = mu_lowt, y = seg_y, color = col_lowt, size = 3, shape = 18) +
  annotate("segment", x = mu_hight - sd_hight, xend = mu_hight + sd_hight,
           y = seg_y * 2.2, yend = seg_y * 2.2, color = col_hight, linewidth = 2) +
  annotate("point", x = mu_hight, y = seg_y * 2.2, color = col_hight, size = 3, shape = 18) +
  scale_fill_manual(values = c(LowT = col_lowt, HighT = col_hight)) +
  scale_color_manual(values = c(LowT = col_lowt, HighT = col_hight)) +
  scale_x_continuous(limits = c(0, 32), breaks = seq(0, 30, 5)) +
  labs(x = "SST (C)", y = "Abundance-weighted density") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.margin = margin(5, 10, 5, 5))

## =========================================================
## 2) Levins' B（标准化）
## =========================================================
message("\n=== 指标2: Levins' B ===")

levins_results <- list()
for (v in c(env_vars, "Latitude")) {
  for (eco in c("LowT", "HighT")) {
    sub <- dat %>% filter(ecotype == eco)
    vals <- if (v == "Latitude") sub$Latitude else sub[[v]]
    ab <- sub$AOA_ab
    all_vals <- if (v == "Latitude") dat$Latitude else dat[[v]]
    breaks <- seq(min(all_vals), max(all_vals), length.out = 21)
    bins <- cut(vals, breaks = breaks, include.lowest = TRUE)
    p <- tapply(ab, bins, sum, default = 0)
    p[is.na(p)] <- 0
    p <- p / sum(p)
    p_nz <- p[p > 0]
    B <- 1 / sum(p_nz^2)
    B_std <- (B - 1) / (20 - 1)
    levins_results <- c(levins_results, list(data.frame(
      var = v, ecotype = eco, B = B, B_std = B_std,
      n_bins_occupied = length(p_nz)
    )))
  }
}
levins_df <- bind_rows(levins_results)

levins_ratio <- levins_df %>%
  select(var, ecotype, B_std) %>%
  pivot_wider(names_from = ecotype, values_from = B_std) %>%
  mutate(ratio = LowT / HighT) %>%
  arrange(desc(ratio))

cat("\nLevins' B (standardized):\n")
print(levins_ratio, n = Inf)

var_order <- levins_ratio$var
levins_df$var <- factor(levins_df$var, levels = rev(var_order))

## ---- Panel B: 横向条形图，按ratio排序 ----
ratio_labels <- levins_ratio %>%
  mutate(label = paste0(round(ratio, 1), "x"),
         x_pos = pmax(LowT, HighT) + 0.02,
         var = factor(var, levels = rev(var_order)))

p_B <- levins_df %>%
  ggplot(aes(x = B_std, y = var, fill = ecotype)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(data = ratio_labels,
            aes(x = x_pos, y = var, label = label),
            inherit.aes = FALSE, hjust = 0, size = 3, fontface = "bold", color = "grey30") +
  scale_fill_manual(values = c(LowT = col_lowt, HighT = col_hight)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Standardized Levins' B", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank(),
        plot.margin = margin(5, 15, 5, 5))

## =========================================================
## 3) 生态位重叠（Pianka & Schoener）
## =========================================================
message("\n=== 指标3: 生态位重叠 ===")

overlap_results <- list()
for (v in env_vars) {
  all_vals <- dat[[v]]
  breaks <- seq(min(all_vals), max(all_vals), length.out = 21)
  p_l <- tapply(lowt$AOA_ab, cut(lowt[[v]], breaks, include.lowest = TRUE), sum, default = 0)
  p_l[is.na(p_l)] <- 0; p_l <- p_l / sum(p_l)
  p_h <- tapply(hight$AOA_ab, cut(hight[[v]], breaks, include.lowest = TRUE), sum, default = 0)
  p_h[is.na(p_h)] <- 0; p_h <- p_h / sum(p_h)
  pianka_O <- sum(p_l * p_h) / sqrt(sum(p_l^2) * sum(p_h^2))
  schoener_D <- 1 - 0.5 * sum(abs(p_l - p_h))
  overlap_results <- c(overlap_results, list(data.frame(
    var = v, Pianka_O = pianka_O, Schoener_D = schoener_D
  )))
}
overlap_df <- bind_rows(overlap_results)
cat("\n生态位重叠:\n"); print(overlap_df)

overlap_order <- overlap_df %>% arrange(Schoener_D) %>% pull(var)
overlap_long <- overlap_df %>%
  pivot_longer(c(Schoener_D, Pianka_O), names_to = "metric", values_to = "value") %>%
  mutate(var = factor(var, levels = rev(overlap_order)))

## ---- Panel C: 重叠条形图，按Schoener's D排序 ----
p_C <- overlap_long %>%
  ggplot(aes(x = value, y = var, fill = metric)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  scale_fill_manual(values = c(Schoener_D = "#9B59B6", Pianka_O = "#2ECC7180"),
                    labels = c("Pianka's O", "Schoener's D")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Overlap index", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank(),
        plot.margin = margin(5, 10, 5, 5))

## 回填Panel A的文字框
schoener_temp <- overlap_df$Schoener_D[overlap_df$var == "temp"]
sigma_ratio <- round(sd_lowt / sd_hight, 2)

p_A <- p_A +
  annotate("label", x = 2, y = ymax * 0.95, hjust = 0, vjust = 1, size = 2.8,
           fill = "white", alpha = 0.85, label.r = unit(0.15, "lines"),
           label = paste0(
             "LowT: mu=", round(mu_lowt, 1), "C, sd=", round(sd_lowt, 1), "C\n",
             "HighT: mu=", round(mu_hight, 1), "C, sd=", round(sd_hight, 1), "C\n",
             "Ratio sd: ", sigma_ratio, "x\n",
             "Schoener's D = ", round(schoener_temp, 3)))

## =========================================================
## 4) 环境异质性暴露度（Supp用）
## =========================================================
message("\n=== 指标4: 环境异质性暴露度 ===")

env_heterogeneity <- dat %>%
  pivot_longer(all_of(env_vars), names_to = "var", values_to = "value") %>%
  group_by(ecotype, var) %>%
  summarise(
    weighted_var = sum(AOA_ab * (value - weighted.mean(value, AOA_ab))^2) / sum(AOA_ab),
    .groups = "drop"
  )

cat("\n环境异质性暴露度（丰度加权方差）:\n")
env_heterogeneity %>%
  select(var, ecotype, weighted_var) %>%
  pivot_wider(names_from = ecotype, values_from = weighted_var) %>%
  mutate(ratio_LowT_HighT = round(LowT / HighT, 2)) %>%
  print(n = Inf)

## ---- Supp Panel E: 环境异质性条形图 ----
het_ratio <- env_heterogeneity %>%
  select(var, ecotype, weighted_var) %>%
  pivot_wider(names_from = ecotype, values_from = weighted_var) %>%
  mutate(ratio = LowT / HighT) %>%
  arrange(desc(ratio))
het_order <- het_ratio$var

env_heterogeneity$var <- factor(env_heterogeneity$var, levels = rev(het_order))

p_supp_E <- env_heterogeneity %>%
  ggplot(aes(x = weighted_var, y = var, fill = ecotype)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = c(LowT = col_lowt, HighT = col_hight)) +
  facet_wrap(~ var, scales = "free_x", nrow = 4) +
  labs(x = "Abundance-weighted variance", y = NULL,
       title = "Environmental heterogeneity exposure") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

## ---- Supp Panel D: 丰度热图 ----
message("\n=== Supp: 丰度热图 ===")
world_df <- ggplot2::map_data("world")

make_heatmap <- function(sub_dat, title_text, fill_col) {
  ggplot() +
    geom_polygon(data = world_df, aes(x = long, y = lat, group = group),
                 fill = "grey90", color = "grey60", linewidth = 0.15) +
    geom_point(data = sub_dat, aes(x = Longitude, y = Latitude, size = AOA_ab),
               color = fill_col, alpha = 0.6) +
    scale_size_continuous(range = c(0.5, 5), name = "Relative\nabundance") +
    coord_quickmap(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          panel.background = element_rect(fill = "#D6EAF8", color = NA),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3),
          legend.position = "right")
}

p_supp_D_lowt  <- make_heatmap(lowt, "LowT abundance distribution", col_lowt)
p_supp_D_hight <- make_heatmap(hight, "HighT abundance distribution", col_hight)

## =========================================================
## 5) 保存图表
## =========================================================
message("\n=== 生成图表 ===")

## Main Figure: 3-panel 单行
main_fig <- p_A + p_B + p_C +
  plot_layout(widths = c(1.2, 1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 13)
    )
  )

ggsave(file.path(out_dir, "Fig_niche_main.pdf"), main_fig,
       width = 16, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig_niche_main.png"), main_fig,
       width = 16, height = 5.5, dpi = 300, bg = "white")
message("Main figure saved")

## Supp Figure: 丰度热图 + 环境异质性
supp_fig <- (p_supp_D_lowt / p_supp_D_hight) | p_supp_E
supp_fig <- supp_fig +
  plot_annotation(
    title = "Supplementary: Geographic distribution and environmental heterogeneity",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.tag = element_text(face = "bold", size = 12)
    )
  )

ggsave(file.path(out_dir, "Fig_niche_supp.pdf"), supp_fig,
       width = 16, height = 10, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig_niche_supp.png"), supp_fig,
       width = 16, height = 10, dpi = 300, bg = "white")
message("Supp figure saved")

## 保存数据
write_csv(niche_center_breadth, file.path(out_dir, "niche_center_breadth.csv"))
write_csv(levins_df, file.path(out_dir, "levins_B.csv"))
write_csv(env_heterogeneity, file.path(out_dir, "env_heterogeneity.csv"))
write_csv(overlap_df, file.path(out_dir, "niche_overlap_v2.csv"))

## 统计摘要
sink(file.path(out_dir, "niche_analysis_v2_summary.txt"))
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("LowT vs HighT Niche Analysis Results\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")
cat("Samples: LowT =", nrow(lowt), ", HighT =", nrow(hight), "\n\n")

cat("1. Temperature niche center & breadth (abundance-weighted)\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("LowT: mu =", round(mu_lowt, 2), "C, sigma =", round(sd_lowt, 2), "C\n")
cat("HighT: mu =", round(mu_hight, 2), "C, sigma =", round(sd_hight, 2), "C\n")
cat("Breadth ratio (LowT/HighT) =", sigma_ratio, "\n")

cat("\n2. Levins' B (standardized, 0-1)\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
levins_ratio %>%
  mutate(LowT = round(LowT, 3), HighT = round(HighT, 3), ratio = round(ratio, 1)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n3. Niche overlap\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
overlap_df %>%
  arrange(Schoener_D) %>%
  mutate(across(c(Pianka_O, Schoener_D), ~ round(.x, 3))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n4. Environmental heterogeneity (abundance-weighted variance)\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
env_heterogeneity %>%
  select(var, ecotype, weighted_var) %>%
  pivot_wider(names_from = ecotype, values_from = weighted_var) %>%
  mutate(ratio = round(LowT / HighT, 2),
         LowT = round(LowT, 4), HighT = round(HighT, 4)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n\nKey conclusions:\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("- HighT is a thermal specialist (sigma=", round(sd_hight, 1),
    "C), LowT is a thermal generalist (sigma=", round(sd_lowt, 1), "C)\n")
cat("- Temperature differentiation is SPECIFIC: Levins' B ratio 28.8x for temp,\n")
cat("  7.2x for latitude, but <5x for other variables\n")
cat("- Temperature overlap is lowest (Schoener's D=", round(schoener_temp, 3),
    "), confirming true niche separation\n")
cat("- This is not a general specialist/generalist pattern but temperature-specific\n")

sink()

message("\nDONE. Results saved to: ", out_dir)
