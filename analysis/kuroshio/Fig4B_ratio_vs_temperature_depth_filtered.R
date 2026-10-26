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
# Fig4B: Ratio vs Temperature scatter plot
# 修正版：添加深度过滤（0-200m）
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
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
  filter(is.finite(Temperature), is.finite(ratio))

# 读取 water mass 标签
wm_csv <- file.path(workdir, "output_tables", "watermass_labels.csv")
if (file.exists(wm_csv)) {
  wm_labels <- read_csv(wm_csv, show_col_types = FALSE)
  dat <- dat %>%
    left_join(wm_labels %>% select(Station, WaterMass), by = "Station")
}

# ---- 深度过滤：只保留 0-200m ----
cat("\n=== Depth filtering ===\n")
cat("Total samples before filtering:", nrow(dat), "\n")
cat("Depth range:", min(dat$Depth_m, na.rm = TRUE), "-", max(dat$Depth_m, na.rm = TRUE), "m\n")

dat <- dat %>%
  filter(Depth_m <= 200)

cat("Samples after filtering (0-200m):", nrow(dat), "\n")
cat("Depth range after filtering:", min(dat$Depth_m, na.rm = TRUE), "-", max(dat$Depth_m, na.rm = TRUE), "m\n\n")

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

# ---- 2. 计算相关性 ----
cor_result <- cor.test(dat$Temperature, dat$ratio, method = "spearman")
cor_label <- sprintf("Spearman ρ = %.3f\np = %.3f", cor_result$estimate, cor_result$p.value)

# ---- 3. Style settings ----
wm_colors <- c(
  "Oyashio"    = "#4C78A8",
  "Transition" = "grey70",
  "Kuroshio"   = "#F58518"
)

wm_shapes <- c(
  "Oyashio"    = 18,
  "Transition" = 15,
  "Kuroshio"   = 16
)

# ---- 4. 绘图 ----
p <- ggplot(dat, aes(x = Temperature, y = ratio)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey60") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80",
              linewidth = 1, alpha = 0.3) +
  geom_point(aes(color = WaterMass, shape = WaterMass), size = 6, alpha = 0.85) +
  scale_color_manual(values = wm_colors, name = "Water mass") +
  scale_shape_manual(values = wm_shapes, name = "Water mass") +
  annotate("text", x = min(dat$Temperature) + 1, y = max(dat$ratio) - 0.2,
           label = cor_label, hjust = 0, size = 4, color = "grey20") +
  labs(
    title = "Ratio vs Temperature (depth-filtered: 0-200m)",
    x = "Temperature (°C)",
    y = expression(log[10](HighT/LowT))
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 11),
    legend.position = "right",
    plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.margin  = margin(10, 10, 10, 10)
  )

# ---- 5. 保存 ----
ggsave(file.path(workdir, "Fig4B_ratio_vs_temperature_simplified.pdf"), p,
       width = 6.9, height = 6)
ggsave(file.path(workdir, "Fig4B_ratio_vs_temperature_simplified.png"), p,
       width = 8, height = 6, dpi = 400)

cat("\n✓ Fig4B (depth-filtered, simplified classification) saved to:", workdir, "\n")
cat("  - Fig4B_ratio_vs_temperature_simplified.pdf\n")
cat("  - Fig4B_ratio_vs_temperature_simplified.png\n\n")

# ---- 6. 统计摘要 ----
cat("Statistical summary:\n")
cat("Spearman correlation: ρ =", round(cor_result$estimate, 3), ", p =",
    formatC(cor_result$p.value, format = "f", digits = 4), "\n")
cat("\nSample counts by water mass:\n")
print(table(dat$WaterMass))
cat("\nTemperature range:", min(dat$Temperature), "-", max(dat$Temperature), "°C\n")
cat("Ratio range:", round(min(dat$ratio), 2), "-", round(max(dat$ratio), 2), "\n")
