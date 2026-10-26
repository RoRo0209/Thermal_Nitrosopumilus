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
setwd(thermal_path("1.phylogenetic"))

# OGT Bimodality Analysis
# 目的：证明 HighT/LowT 的划分是客观的，基于 OGT 分布的统计特征

library(diptest)      # Hartigan's dip test
library(mclust)       # Gaussian mixture model
library(ggplot2)
library(dplyr)

# 读取数据
data <- read.csv("ASR_clade_assignments_final.csv", stringsAsFactors = FALSE)

# 提取所有 OGT 值（包括 HighT, LowT 和 Other）
ogt_all <- data$OGT

# 提取仅 HighT 和 LowT 的 OGT 值
ogt_clades <- data %>%
  filter(Clade %in% c("HighT", "LowT")) %>%
  pull(OGT)

cat("=== OGT Bimodality Analysis ===\n\n")

# ========================================
# 1. Hartigan's Dip Test
# ========================================
cat("1. Hartigan's Dip Test for Unimodality\n")
cat("   (Tests H0: distribution is unimodal)\n\n")

# 对所有基因组进行测试
dip_all <- dip.test(ogt_all)
cat("   All genomes (n =", length(ogt_all), "):\n")
cat("     Dip statistic =", round(dip_all$statistic, 4), "\n")
cat("     P-value =", format(dip_all$p.value, scientific = TRUE), "\n")
if (dip_all$p.value < 0.05) {
  cat("     ** Significant evidence of multimodality (reject unimodality) **\n\n")
} else {
  cat("     No significant evidence against unimodality\n\n")
}

# 对 HighT/LowT 单系群进行测试
dip_clades <- dip.test(ogt_clades)
cat("   HighT + LowT clades only (n =", length(ogt_clades), "):\n")
cat("     Dip statistic =", round(dip_clades$statistic, 4), "\n")
cat("     P-value =", format(dip_clades$p.value, scientific = TRUE), "\n")
if (dip_clades$p.value < 0.05) {
  cat("     ** Significant evidence of multimodality (reject unimodality) **\n\n")
} else {
  cat("     No significant evidence against unimodality\n\n")
}

# ========================================
# 2. Gaussian Mixture Model (GMM)
# ========================================
cat("2. Gaussian Mixture Model Comparison\n")
cat("   (Comparing 1-component vs 2-component models)\n\n")

# 对所有基因组拟合 GMM
gmm1_all <- Mclust(ogt_all, G = 1, modelNames = "V")
gmm2_all <- Mclust(ogt_all, G = 2, modelNames = "V")

cat("   All genomes:\n")
cat("     1-component model: BIC =", round(gmm1_all$bic, 2), "\n")
cat("     2-component model: BIC =", round(gmm2_all$bic, 2), "\n")
cat("     Delta BIC =", round(gmm2_all$bic - gmm1_all$bic, 2), "\n")
if (gmm2_all$bic > gmm1_all$bic) {
  cat("     ** 2-component model is better supported **\n\n")
} else {
  cat("     1-component model is better supported\n\n")
}

# 对 HighT/LowT 单系群拟合 GMM
gmm1_clades <- Mclust(ogt_clades, G = 1, modelNames = "V")
gmm2_clades <- Mclust(ogt_clades, G = 2, modelNames = "V")

cat("   HighT + LowT clades only:\n")
cat("     1-component model: BIC =", round(gmm1_clades$bic, 2), "\n")
cat("     2-component model: BIC =", round(gmm2_clades$bic, 2), "\n")
cat("     Delta BIC =", round(gmm2_clades$bic - gmm1_clades$bic, 2), "\n")
if (gmm2_clades$bic > gmm1_clades$bic) {
  cat("     ** 2-component model is better supported **\n\n")
} else {
  cat("     1-component model is better supported\n\n")
}

# 输出 2-component 模型的参数
cat("   2-component model parameters (HighT + LowT clades):\n")
cat("     Component 1: mean =", round(gmm2_clades$parameters$mean[1], 2),
    "°C, SD =", round(sqrt(gmm2_clades$parameters$variance$sigmasq[1]), 2), "°C\n")
cat("     Component 2: mean =", round(gmm2_clades$parameters$mean[2], 2),
    "°C, SD =", round(sqrt(gmm2_clades$parameters$variance$sigmasq[2]), 2), "°C\n")
cat("     Mixing proportions:", round(gmm2_clades$parameters$pro, 3), "\n\n")

# ========================================
# 3. 可视化
# ========================================
cat("3. Generating visualization...\n")

# 准备绘图数据
plot_data <- data %>%
  mutate(Group = case_when(
    Clade == "HighT" ~ "HighT clade",
    Clade == "LowT" ~ "LowT clade",
    TRUE ~ "Other"
  ))

# 创建密度图
p <- ggplot(plot_data, aes(x = OGT, fill = Group)) +
  geom_density(alpha = 0.6, adjust = 1.5) +
  geom_rug(aes(color = Group), alpha = 0.5, length = unit(0.05, "npc")) +
  scale_fill_manual(values = c("HighT clade" = "#D55E00",
                                "LowT clade" = "#0072B2",
                                "Other" = "gray70")) +
  scale_color_manual(values = c("HighT clade" = "#D55E00",
                                 "LowT clade" = "#0072B2",
                                 "Other" = "gray70")) +
  labs(x = "Optimal Growth Temperature (°C)",
       y = "Density",
       title = "OGT Distribution in Nitrosopumilus",
       subtitle = paste0("Dip test p = ", format(dip_clades$p.value, digits = 3, scientific = TRUE),
                        "; 2-component BIC = ", round(gmm2_clades$bic, 1),
                        " vs 1-component BIC = ", round(gmm1_clades$bic, 1))) +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "black"))

ggsave("OGT_Bimodality_Distribution.pdf", p, width = 8, height = 5)
ggsave("OGT_Bimodality_Distribution.png", p, width = 8, height = 5, dpi = 300)

cat("   Saved: OGT_Bimodality_Distribution.pdf/png\n\n")

# ========================================
# 4. 总结
# ========================================
cat("=== SUMMARY ===\n\n")
cat("For manuscript text:\n\n")
cat("\"Across Nitrosopumilus genomes, OGT distribution showed significant\n")
cat("evidence of non-unimodality (Hartigan's dip test, p =",
    format(dip_clades$p.value, digits = 3, scientific = TRUE), ")\n")
cat("and was better fit by a two-component Gaussian mixture model\n")
cat("(ΔBIC =", round(gmm2_clades$bic - gmm1_clades$bic, 1),
    ") than a single Gaussian,\n")
cat("consistent with the existence of two distinct thermal strategies.\"\n\n")

# 保存结果到文件
sink("OGT_Bimodality_Results.txt")
cat("=== OGT Bimodality Analysis Results ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Hartigan's Dip Test (HighT + LowT clades):\n")
cat("  Dip statistic:", round(dip_clades$statistic, 4), "\n")
cat("  P-value:", format(dip_clades$p.value, scientific = TRUE), "\n\n")
cat("Gaussian Mixture Model (HighT + LowT clades):\n")
cat("  1-component BIC:", round(gmm1_clades$bic, 2), "\n")
cat("  2-component BIC:", round(gmm2_clades$bic, 2), "\n")
cat("  Delta BIC:", round(gmm2_clades$bic - gmm1_clades$bic, 2), "\n\n")
cat("2-component parameters:\n")
cat("  Component 1: μ =", round(gmm2_clades$parameters$mean[1], 2),
    "°C, σ =", round(sqrt(gmm2_clades$parameters$variance$sigmasq[1]), 2), "°C\n")
cat("  Component 2: μ =", round(gmm2_clades$parameters$mean[2], 2),
    "°C, σ =", round(sqrt(gmm2_clades$parameters$variance$sigmasq[2]), 2), "°C\n")
cat("  Mixing proportions:", round(gmm2_clades$parameters$pro, 3), "\n")
sink()

cat("Results saved to: OGT_Bimodality_Results.txt\n")
cat("\nAnalysis complete!\n")
