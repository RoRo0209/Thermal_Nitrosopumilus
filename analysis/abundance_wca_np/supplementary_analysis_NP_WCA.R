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
setwd(file.path(.release_root, "abundance_wca_np"))

## =========================================================
## supplementary_analysis_NP_WCA.R
## 补充分析：按优先.txt要求（仅NP vs WCA）
##   1. 温度对NP/WCA比例的偏相关分析（控制其他环境变量）
##   2. NP/WCA对温度的响应曲线（GAM + 置信区间）
##   3. 温度×营养盐交互项分析
## =========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)      # GAM
  library(ppcor)     # 偏相关
  library(patchwork)
})

# 确保使用dplyr的select
select <- dplyr::select

## ---- 路径 ----
base_dir <- file.path(.release_root, "abundance_wca_np")
np_wca_file <- file.path(base_dir, "NP_WCA_final.xlsx")
geo_file    <- file.path(base_dir, "genome_lon_lat.xlsx")

biooracle_grid <- .release_input("THERMAL_BASELINE_GRID")

out_dir <- file.path(base_dir, "NP_WCA_global_ratio")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grid_res_deg <- 0.25

## ---- 辅助函数 ----
align_lonlat_to_grid <- function(df, res_deg,
                                 xlim = c(-180, 180),
                                 ylim = c(-90,  90)) {
  df %>%
    mutate(
      Longitude = floor((.data[["Longitude"]] - xlim[1]) / res_deg) * res_deg + xlim[1] + res_deg / 2,
      Latitude  = floor((.data[["Latitude"]]  - ylim[1]) / res_deg) * res_deg + ylim[1] + res_deg / 2
    )
}

## =========================================================
## 1) 读取数据并提取环境变量
## =========================================================
message("=== 读取NP/WCA数据 ===")

## NP/WCA
np_wca <- read_excel(np_wca_file) %>%
  rename(sample = 1, NP_frac = 2, WCA_frac = 3) %>%
  filter(is.finite(NP_frac), is.finite(WCA_frac)) %>%
  mutate(
    total = NP_frac + WCA_frac,
    NP_frac = NP_frac / total,
    WCA_frac = WCA_frac / total
  ) %>%
  select(-total)

## 地理坐标
geo <- read_excel(geo_file) %>%
  rename(sample = 1, Longitude = 2, Latitude = 3) %>%
  mutate(
    Longitude = suppressWarnings(as.numeric(Longitude)),
    Latitude  = suppressWarnings(as.numeric(Latitude))
  ) %>%
  filter(is.finite(Longitude), is.finite(Latitude))

## 合并
dat_np <- inner_join(np_wca, geo, by = "sample")
message("NP/WCA样本数: ", nrow(dat_np))

## 提取环境变量
message("读取Bio-ORACLE环境变量...")
base_grid <- read_tsv(biooracle_grid, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    Longitude = suppressWarnings(as.numeric(Longitude)),
    Latitude  = suppressWarnings(as.numeric(Latitude))
  ) %>%
  filter(is.finite(Longitude), is.finite(Latitude))

pred_vars <- c("temp", "sal", "chl", "si", "no3", "po4", "o2")
pred_vars <- intersect(pred_vars, names(base_grid))

base_grid_aligned <- base_grid %>%
  align_lonlat_to_grid(grid_res_deg) %>%
  group_by(Longitude, Latitude) %>%
  summarise(across(all_of(pred_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

## 为样本提取环境变量
dat_np_aligned <- dat_np %>% align_lonlat_to_grid(grid_res_deg)
dat_np_env <- dat_np_aligned %>%
  left_join(base_grid_aligned, by = c("Longitude", "Latitude")) %>%
  drop_na(all_of(pred_vars))

message("NP/WCA有效样本（含环境变量）: ", nrow(dat_np_env))

## =========================================================
## 2) 分析1: 温度对NP/WCA比例的偏相关（控制其他变量）
## =========================================================
message("\n=== 分析1: 偏相关分析 ===")

## 原始Spearman相关
cor_raw_np <- cor.test(dat_np_env$temp, dat_np_env$NP_frac, method = "spearman")
cat("\nNP/WCA vs 温度（原始Spearman相关）:\n")
cat("  ρ = ", round(cor_raw_np$estimate, 3), ", p = ", format.pval(cor_raw_np$p.value, digits = 3), "\n")

## 偏Spearman相关：控制 sal, no3, po4, si, o2
control_vars <- c("sal", "no3", "po4", "si", "o2")
control_vars <- intersect(control_vars, names(dat_np_env))

if (length(control_vars) >= 2) {
  pcor_np <- pcor.test(
    dat_np_env$NP_frac,
    dat_np_env$temp,
    dat_np_env[, control_vars],
    method = "spearman"
  )
  cat("\nNP/WCA vs 温度（偏Spearman相关，控制", paste(control_vars, collapse=", "), "）:\n")
  cat("  ρ_partial = ", round(pcor_np$estimate, 3), ", p = ", format.pval(pcor_np$p.value, digits = 3), "\n")

  cat("\n解释：\n")
  if (sign(cor_raw_np$estimate) != sign(pcor_np$estimate)) {
    cat("  ⚠️  控制其他变量后，温度效应发生反转！\n")
    cat("  原始负相关可能是混杂效应（低温水体通常高纬度、高营养、高DO）\n")
  } else if (abs(pcor_np$estimate) < abs(cor_raw_np$estimate) * 0.5) {
    cat("  ⚠️  控制其他变量后，温度效应大幅减弱\n")
    cat("  说明温度的表观效应部分来自与其他因子的共线性\n")
  } else {
    cat("  ✓  温度效应在控制其他变量后仍然稳健\n")
  }
}

## =========================================================
## 3) 分析2: GAM响应曲线（温度 → NP/WCA）
## =========================================================
message("\n=== 分析2: GAM温度响应曲线 ===")

## 单变量GAM: NP/WCA ~ temp
gam_np_temp <- gam(NP_frac ~ s(temp, k = 8), data = dat_np_env, method = "REML")
cat("\n单变量GAM: NP/WCA ~ s(temp)\n")
cat("  R² = ", round(summary(gam_np_temp)$r.sq, 3), "\n")
cat("  Deviance explained = ", round(summary(gam_np_temp)$dev.expl * 100, 1), "%\n")
cat("  温度项 p-value: ", format.pval(summary(gam_np_temp)$s.table[1, "p-value"]), "\n")

## 预测 + 置信区间
temp_seq <- seq(min(dat_np_env$temp, na.rm = TRUE),
                max(dat_np_env$temp, na.rm = TRUE),
                length.out = 200)
pred_np <- predict(gam_np_temp, newdata = data.frame(temp = temp_seq), se.fit = TRUE)
pred_np_df <- data.frame(
  temp = temp_seq,
  fit = pred_np$fit,
  se = pred_np$se.fit,
  lower = pred_np$fit - 1.96 * pred_np$se.fit,
  upper = pred_np$fit + 1.96 * pred_np$se.fit
)

## 可视化
p_gam_np <- ggplot() +
  geom_point(data = dat_np_env, aes(x = temp, y = NP_frac),
             alpha = 0.3, size = 1.5, color = "grey50") +
  geom_ribbon(data = pred_np_df, aes(x = temp, ymin = lower, ymax = upper),
              fill = "#E8A07A", alpha = 0.3) +
  geom_line(data = pred_np_df, aes(x = temp, y = fit),
            color = "#E8A07A", linewidth = 1.5) +
  labs(x = "Temperature (°C)", y = "NP fraction",
       title = "NP vs WCA: Temperature response (GAM with 95% CI)",
       subtitle = paste0("R² = ", round(summary(gam_np_temp)$r.sq, 3),
                        ", Deviance explained = ", round(summary(gam_np_temp)$dev.expl * 100, 1), "%")) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10))

## =========================================================
## 4) 分析3: 温度×营养盐交互项
## =========================================================
message("\n=== 分析3: 温度×营养盐交互项 ===")

## NP/WCA ~ temp × no3 + 其他变量
if ("no3" %in% names(dat_np_env)) {
  lm_np_interact <- lm(NP_frac ~ temp * no3 + sal + po4 + si + o2, data = dat_np_env)
  cat("\n线性模型: NP/WCA ~ temp × NO3 + 其他变量\n")
  cat("  R² = ", round(summary(lm_np_interact)$r.squared, 3), "\n")
  cat("  Adjusted R² = ", round(summary(lm_np_interact)$adj.r.squared, 3), "\n")

  ## 提取交互项系数
  coef_interact <- coef(summary(lm_np_interact))
  if ("temp:no3" %in% rownames(coef_interact)) {
    cat("\n交互项 temp:no3:\n")
    cat("  β = ", round(coef_interact["temp:no3", "Estimate"], 5),
        ", p = ", format.pval(coef_interact["temp:no3", "Pr(>|t|)"], digits = 3), "\n")

    if (coef_interact["temp:no3", "Pr(>|t|)"] < 0.05) {
      cat("  ✓ 交互项显著：温度和NO3对NP/WCA有协同效应\n")
    } else {
      cat("  × 交互项不显著：温度和NO3效应相互独立\n")
    }
  }

  ## 可视化交互效应
  # 分3个NO3水平展示温度效应
  no3_levels <- quantile(dat_np_env$no3, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
  pred_interact <- expand.grid(
    temp = temp_seq,
    no3 = no3_levels,
    sal = mean(dat_np_env$sal, na.rm = TRUE),
    po4 = mean(dat_np_env$po4, na.rm = TRUE),
    si = mean(dat_np_env$si, na.rm = TRUE),
    o2 = mean(dat_np_env$o2, na.rm = TRUE)
  )
  pred_interact$NP_frac <- predict(lm_np_interact, newdata = pred_interact)
  pred_interact$NO3_level <- factor(pred_interact$no3,
                                     labels = c("Low NO3 (10%ile)", "Medium NO3 (50%ile)", "High NO3 (90%ile)"))

  p_interact <- ggplot(pred_interact, aes(x = temp, y = NP_frac, color = NO3_level)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("#3498DB", "#F39C12", "#E74C3C")) +
    labs(x = "Temperature (°C)", y = "Predicted NP fraction",
         title = "Temperature × NO3 interaction effect on NP/WCA",
         color = NULL) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "top")
}

## =========================================================
## 5) 多变量GAM（控制其他变量后的温度效应）
## =========================================================
message("\n=== 分析4: 多变量GAM（控制其他环境变量）===")

## NP/WCA ~ s(temp) + s(sal) + s(no3) + ...
gam_np_multi <- gam(
  NP_frac ~ s(temp, k = 6) + s(sal, k = 5) + s(no3, k = 5) + s(po4, k = 5) + s(si, k = 5) + s(o2, k = 5),
  data = dat_np_env,
  method = "REML"
)
cat("\n多变量GAM: NP/WCA ~ s(temp) + s(其他变量)\n")
cat("  R² = ", round(summary(gam_np_multi)$r.sq, 3), "\n")
cat("  Deviance explained = ", round(summary(gam_np_multi)$dev.expl * 100, 1), "%\n")
cat("\n各变量显著性:\n")
print(summary(gam_np_multi)$s.table)

## 提取温度的偏效应曲线
cat("\n温度的偏效应（控制其他变量）:\n")
cat("  edf = ", round(summary(gam_np_multi)$s.table[1, "edf"], 2), "\n")
cat("  p-value = ", format.pval(summary(gam_np_multi)$s.table[1, "p-value"]), "\n")

## 可视化温度偏效应
temp_partial <- plot(gam_np_multi, select = 1, seWithMean = TRUE)
temp_partial_df <- data.frame(
  temp = temp_partial[[1]]$x,
  fit = temp_partial[[1]]$fit,
  se = temp_partial[[1]]$se
)
temp_partial_df$lower <- temp_partial_df$fit - 1.96 * temp_partial_df$se
temp_partial_df$upper <- temp_partial_df$fit + 1.96 * temp_partial_df$se

p_partial <- ggplot(temp_partial_df, aes(x = temp, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#3498DB", alpha = 0.3) +
  geom_line(color = "#3498DB", linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "Temperature (°C)", y = "Partial effect on NP fraction",
       title = "Temperature partial effect (controlling other variables)",
       subtitle = paste0("From multivariate GAM, R² = ", round(summary(gam_np_multi)$r.sq, 3))) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10))

## =========================================================
## 6) 输出图表
## =========================================================
message("\n=== 生成图表 ===")

## 组合图
if (exists("p_interact")) {
  combined <- (p_gam_np | p_partial) / (p_interact | plot_spacer()) +
    plot_annotation(
      title = "NP vs WCA: Temperature effects and interactions",
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.tag = element_text(face = "bold", size = 12)
      )
    )
} else {
  combined <- (p_gam_np | p_partial) +
    plot_annotation(
      title = "NP vs WCA: Temperature effects",
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.tag = element_text(face = "bold", size = 12)
      )
    )
}

out_pdf <- file.path(out_dir, "supplementary_NP_WCA_temperature_analysis.pdf")
ggsave(out_pdf, combined, width = 14, height = 10, dpi = 300, bg = "white")
message("PDF保存: ", out_pdf)

## 保存数据
write_tsv(pred_np_df, file.path(out_dir, "GAM_NP_WCA_temperature_response.tsv"))
write_tsv(temp_partial_df, file.path(out_dir, "GAM_NP_WCA_temperature_partial_effect.tsv"))

## 保存统计结果摘要
sink(file.path(out_dir, "supplementary_NP_WCA_statistics.txt"))
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("NP vs WCA 补充统计分析结果\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("样本数: ", nrow(dat_np_env), "\n")
cat("温度范围: ", round(min(dat_np_env$temp), 1), " ~ ", round(max(dat_np_env$temp), 1), "°C\n")
cat("NP fraction 范围: ", round(min(dat_np_env$NP_frac), 3), " ~ ", round(max(dat_np_env$NP_frac), 3), "\n\n")

cat("1. 偏相关分析\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("原始 Spearman ρ = ", round(cor_raw_np$estimate, 3), ", p = ", format.pval(cor_raw_np$p.value), "\n")
if (exists("pcor_np")) {
  cat("偏 Spearman ρ = ", round(pcor_np$estimate, 3), ", p = ", format.pval(pcor_np$p.value), "\n")
  cat("（控制变量: ", paste(control_vars, collapse = ", "), "）\n")
}

cat("\n2. 单变量GAM\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("R² = ", round(summary(gam_np_temp)$r.sq, 3), "\n")
cat("Deviance explained = ", round(summary(gam_np_temp)$dev.expl * 100, 1), "%\n")
cat("温度项 p-value = ", format.pval(summary(gam_np_temp)$s.table[1, "p-value"]), "\n")

cat("\n3. 多变量GAM\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat("R² = ", round(summary(gam_np_multi)$r.sq, 3), "\n")
cat("Deviance explained = ", round(summary(gam_np_multi)$dev.expl * 100, 1), "%\n")
cat("\n各变量显著性:\n")
print(summary(gam_np_multi)$s.table)

if (exists("lm_np_interact")) {
  cat("\n4. 温度×NO3交互项\n")
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat("模型 R² = ", round(summary(lm_np_interact)$r.squared, 3), "\n")
  if ("temp:no3" %in% rownames(coef_interact)) {
    cat("交互项 temp:no3 β = ", round(coef_interact["temp:no3", "Estimate"], 5),
        ", p = ", format.pval(coef_interact["temp:no3", "Pr(>|t|)"]), "\n")
  }
}

cat("\n" |> rep(2) |> paste(collapse = ""), "\n")
cat("关键结论:\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
if (exists("pcor_np") && sign(cor_raw_np$estimate) != sign(pcor_np$estimate)) {
  cat("• 温度的原始效应是混杂的，控制其他变量后效应反转\n")
  cat("• 溶解氧(O2)是最强预测因子，温度效应较弱\n")
} else {
  cat("• 温度对NP/WCA有独立效应，但解释力有限\n")
}
cat("• 多变量模型可解释", round(summary(gam_np_multi)$dev.expl * 100, 1), "%的变异\n")

sink()

message("\nDONE. 所有结果保存在: ", out_dir)
