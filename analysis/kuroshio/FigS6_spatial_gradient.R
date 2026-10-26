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
setwd(file.path(.release_root, "kuroshio"))

# ============================================================================
# Supplementary Fig. S6: Spatial gradient support
# Panel a: T-S diagram (water mass classification)
# Panel b: End-member comparison (Kuroshio vs Oyashio)
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ---- 0. Paths ----
workdir <- getwd()
in_csv  <- file.path(workdir, "output_tables", "meta_with_env_and_sum.csv")

# ---- 1. Read & prepare data ----
dat <- read_csv(in_csv, show_col_types = FALSE) %>%
  mutate(
    # 计算 ratio = log10((HighT + 1e-6) / (LowT + 1e-6))
    ratio = log10((HighT + 1e-6) / (LowT + 1e-6))
  ) %>%
  filter(is.finite(Temperature), is.finite(Salinity), is.finite(ratio))

# 方案1：只保留表层-次表层（0-200m），避免混淆水团类型和深度效应
# 深层水（>200m）是不同的水团类型，不应与表层黑潮/亲潮比较
cat("\n=== Depth filtering ===\n")
cat("Total samples before filtering:", nrow(dat), "\n")
cat("Depth range:", min(dat$Depth_m, na.rm = TRUE), "-", max(dat$Depth_m, na.rm = TRUE), "m\n")

dat <- dat %>%
  filter(Depth_m <= 200)

cat("Samples after filtering (0-200m):", nrow(dat), "\n")
cat("Depth range after filtering:", min(dat$Depth_m, na.rm = TRUE), "-", max(dat$Depth_m, na.rm = TRUE), "m\n\n")

# 读取 water mass 标签
wm_csv <- file.path(workdir, "output_tables", "watermass_labels.csv")
if (file.exists(wm_csv)) {
  wm_labels <- read_csv(wm_csv, show_col_types = FALSE)
  dat <- dat %>%
    left_join(wm_labels %>% select(Station, WaterMass), by = "Station")
} else {
  # 如果没有预先分类的文件，使用简单的阈值分类
  dat <- dat %>%
    mutate(
      WaterMass = case_when(
        Temperature >= 20 & Salinity >= 34.4 ~ "Kuroshio",
        Temperature >= 15 & Temperature < 20 & Salinity >= 34.4 ~ "Transition",
        Temperature < 10 ~ "Oyashio",
        Temperature < 15 & Salinity < 34.0 ~ "Oyashio",
        TRUE ~ "Transition"
      )
    )
}

# 简化分类：合并 Kuroshio 和 KE_Mixed 为 "Kuroshio"
cat("=== Simplifying water mass classification ===\n")
cat("Original classification:\n")
print(table(dat$WaterMass))

dat <- dat %>%
  mutate(
    WaterMass_simple = case_when(
      WaterMass %in% c("Kuroshio", "KE_Mixed") ~ "Kuroshio",
      WaterMass == "Oyashio" ~ "Oyashio",
      TRUE ~ "Transition"
    ),
    WaterMass = factor(WaterMass_simple, levels = c("Oyashio", "Transition", "Kuroshio"))
  ) %>%
  select(-WaterMass_simple)

cat("\nSimplified classification (Kuroshio + KE_Mixed -> Kuroshio):\n")
print(table(dat$WaterMass))
cat("\n")

# ---- 2. Style settings ----
wm_colors <- c(
  "Oyashio"    = "#4C78A8",
  "Transition" = "grey70",
  "Kuroshio"   = "#F58518"
)

base_theme <- theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 11),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    plot.title   = element_text(face = "bold", size = 12),
    plot.margin  = margin(6, 8, 6, 8)
  )

# ---- 3. Panel a: T-S diagram ----
p_ts <- ggplot(dat, aes(x = Salinity, y = Temperature, color = WaterMass)) +
  geom_point(size = 2.6, alpha = 0.85) +
  scale_color_manual(values = wm_colors, name = "Water mass") +
  labs(
    title = "a  T–S diagram (water mass classification)",
    x = "Salinity (PSU)",
    y = "Temperature (°C)"
  ) +
  base_theme

# ---- 4. Panel b: End-member boxplot ----
# 只比较端元：Kuroshio vs Oyashio
dat_endmember <- dat %>%
  filter(WaterMass %in% c("Kuroshio", "Oyashio")) %>%
  mutate(WaterMass = factor(WaterMass, levels = c("Oyashio", "Kuroshio")))

# Wilcoxon rank-sum test (two-sided)
wilcox_result <- wilcox.test(ratio ~ WaterMass, data = dat_endmember, alternative = "two.sided")
p_value <- wilcox_result$p.value

# 格式化 p 值
p_label <- if (p_value < 0.001) {
  "p < 0.001"
} else {
  sprintf("p = %.3f", p_value)
}

# 计算标注位置
y_max <- max(dat_endmember$ratio, na.rm = TRUE)
y_pos <- y_max + 0.3

p_box <- ggplot(dat_endmember, aes(x = WaterMass, y = ratio, fill = WaterMass)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey60") +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.7, alpha = 0.75) +
  geom_jitter(aes(color = WaterMass), width = 0.12, size = 2.1, alpha = 0.6) +
  annotate("text", x = 1.5, y = y_pos, label = p_label, size = 4, color = "grey20") +
  scale_fill_manual(values = wm_colors) +
  scale_color_manual(values = wm_colors) +
  labs(
    title = "b  End-member contrast of ecotype ratio",
    x = "Water mass",
    y = expression(log[10](HighT/LowT))
  ) +
  base_theme +
  theme(legend.position = "none")

# ---- 5. Assemble & save ----
final <- p_ts + p_box + plot_layout(widths = c(1.05, 0.95))

ggsave(file.path(workdir, "Supplementary_Fig_S6_depth_filtered.pdf"), final,
       width = 11, height = 4.2)
ggsave(file.path(workdir, "Supplementary_Fig_S6_depth_filtered.png"), final,
       width = 11, height = 4.2, dpi = 400)

cat("\n✓ Supplementary Fig. S6 (depth-filtered: 0-200m) saved to:", workdir, "\n")
cat("  - Supplementary_Fig_S6_depth_filtered.pdf\n")
cat("  - Supplementary_Fig_S6_depth_filtered.png\n\n")

# ---- 6. Print statistics ----
cat("Sample counts by water mass:\n")
print(table(dat$WaterMass))

cat("\n\nEnd-member comparison (Kuroshio vs Oyashio):\n")
cat("Kuroshio samples:", sum(dat_endmember$WaterMass == "Kuroshio"), "\n")
cat("Oyashio samples:", sum(dat_endmember$WaterMass == "Oyashio"), "\n")
cat("Wilcoxon test p-value:", formatC(p_value, format = "e", digits = 3), "\n")

# 输出各组的统计摘要
cat("\n\nRatio summary by water mass:\n")
summary_stats <- dat_endmember %>%
  group_by(WaterMass) %>%
  summarize(
    n = n(),
    mean = mean(ratio, na.rm = TRUE),
    median = median(ratio, na.rm = TRUE),
    sd = sd(ratio, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)
