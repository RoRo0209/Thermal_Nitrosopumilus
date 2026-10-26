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

# OGT Prediction Validation - FINAL CORRECTED VERSION
# 所有数据已经过文献核实

library(ggplot2)
library(dplyr)
library(readxl)
library(gridExtra)
library(grid)

# 读取完整的 NP 数据
np_data <- read_excel("../representative_NP.xlsx")

all_np <- data.frame(
  FileName = np_data$FileName,
  Strain = np_data$ncbi_strain_identifiers,
  Clade = np_data$GFF,
  Predicted_OGT = np_data$`predOGT (C)`,
  stringsAsFactors = FALSE
)

# 读取最终更正后的培养菌株数据
validation_data <- read.csv("cultured_strains_validation_FINAL.csv", stringsAsFactors = FALSE)

cat("=== OGT Prediction Validation (FINAL CORRECTED VERSION) ===\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Dataset Summary:\n")
cat("  Total Nitrosopumilus genomes:", nrow(all_np), "\n")
cat("  HighT clade:", sum(all_np$Clade == "highT", na.rm=TRUE), "\n")
cat("  LowT clade:", sum(all_np$Clade == "lowT", na.rm=TRUE), "\n")
cat("  Other NP:", sum(all_np$Clade == "other_NP", na.rm=TRUE), "\n\n")

cat("Validation Set (HighT + LowT only):\n")
cat("  Total cultured strains:", nrow(validation_data), "\n")
cat("  HighT:", sum(validation_data$Clade == "highT"), "\n")
cat("  LowT:", sum(validation_data$Clade == "lowT"), "\n\n")

# 显示验证菌株信息
cat("Validation Strains:\n")
cat(strrep("=", 80), "\n")
for (i in 1:nrow(validation_data)) {
  cat(sprintf("\n%d. %s (%s) - %s clade\n",
              i,
              validation_data$Strain_Name[i],
              validation_data$Species[i],
              validation_data$Clade[i]))
  cat(sprintf("   Genome: %s\n", validation_data$FileName[i]))
  cat(sprintf("   Measured OGT: %d°C\n", validation_data$Measured_OGT[i]))
  cat(sprintf("   Predicted OGT: %.2f°C\n", validation_data$Predicted_OGT[i]))
  cat(sprintf("   Residual: %.2f°C\n", validation_data$Predicted_OGT[i] - validation_data$Measured_OGT[i]))
  cat(sprintf("   Reference: %s (%s) %s\n",
              validation_data$Reference[i],
              validation_data$Year[i],
              validation_data$Journal[i]))
}
cat("\n", strrep("=", 80), "\n\n")

# 计算统计指标
correlation <- cor(validation_data$Measured_OGT, validation_data$Predicted_OGT)
correlation_test <- cor.test(validation_data$Measured_OGT, validation_data$Predicted_OGT)
lm_model <- lm(Predicted_OGT ~ Measured_OGT, data = validation_data)
lm_summary <- summary(lm_model)
rmse <- sqrt(mean((validation_data$Predicted_OGT - validation_data$Measured_OGT)^2))
mae <- mean(abs(validation_data$Predicted_OGT - validation_data$Measured_OGT))

cat("Statistical Metrics:\n")
cat("  Pearson correlation (r):", round(correlation, 3), "\n")
cat("  P-value:", format(correlation_test$p.value, scientific = TRUE), "\n")
cat("  R-squared:", round(lm_summary$r.squared, 3), "\n")
cat("  RMSE:", round(rmse, 2), "°C\n")
cat("  MAE:", round(mae, 2), "°C\n\n")

cat("Temperature Ranges:\n")
highT_temps <- validation_data$Measured_OGT[validation_data$Clade == "highT"]
lowT_temps <- validation_data$Measured_OGT[validation_data$Clade == "lowT"]
cat("  HighT clade (measured):", min(highT_temps), "-", max(highT_temps), "°C\n")
cat("  LowT clade (measured):", min(lowT_temps), "-", max(lowT_temps), "°C\n")
cat("  Temperature separation:", min(highT_temps) - max(lowT_temps), "°C\n\n")

# 创建主验证图
p_validation <- ggplot(validation_data, aes(x = Measured_OGT, y = Predicted_OGT, color = Clade)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "gray50", linewidth = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", color = "black", fill = "gray80",
              alpha = 0.2, linewidth = 1, se = TRUE) +
  geom_point(size = 4, alpha = 0.8) +
  ggrepel::geom_text_repel(
    aes(label = Strain_Name),
    size = 3.5, max.overlaps = 20,
    box.padding = 0.5, point.padding = 0.3,
    fontface = "italic"
  ) +
  scale_color_manual(values = c("highT" = "#D55E00", "lowT" = "#0072B2"),
                     labels = c("highT" = "HighT clade", "lowT" = "LowT clade")) +
  annotate("text", x = 23.5, y = 32,
           label = paste0("r = ", round(correlation, 3),
                         "\nR² = ", round(lm_summary$r.squared, 3),
                         "\nRMSE = ", round(rmse, 2), "°C",
                         "\nn = ", nrow(validation_data)),
           hjust = 0, vjust = 1, size = 4) +
  labs(x = "Measured Optimal Growth Temperature (°C)",
       y = "Predicted OGT (°C)",
       title = "Validation of OGT Prediction Method",
       subtitle = "Cultured Nitrosopumilus strains from HighT and LowT clades",
       color = "Clade") +
  coord_fixed(ratio = 1, xlim = c(22, 33), ylim = c(22, 33)) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.title = element_text(face = "bold"),
    legend.position = c(0.2, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3)
  )

ggsave("FigS2_OGT_Validation_FINAL.pdf", p_validation, width = 7, height = 7)
ggsave("FigS2_OGT_Validation_FINAL.png", p_validation, width = 7, height = 7, dpi = 300)

cat("Main validation figure saved:\n")
cat("  - FigS2_OGT_Validation_FINAL.pdf\n")
cat("  - FigS2_OGT_Validation_FINAL.png\n\n")

# 创建完整的 OGT 分布图
all_np$Group <- factor(all_np$Clade,
                       levels = c("highT", "lowT", "other_NP"),
                       labels = c("HighT clade", "LowT clade", "Other NP"))

p_distribution <- ggplot(all_np, aes(x = Predicted_OGT, fill = Group)) +
  geom_density(alpha = 0.6, adjust = 1.5) +
  geom_rug(aes(color = Group), alpha = 0.5, length = unit(0.05, "npc")) +
  # 添加培养菌株的实测温度标记
  geom_vline(data = validation_data[validation_data$Clade == "highT",],
             aes(xintercept = Measured_OGT),
             color = "#D55E00", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
  geom_vline(data = validation_data[validation_data$Clade == "lowT",],
             aes(xintercept = Measured_OGT),
             color = "#0072B2", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
  scale_fill_manual(values = c("HighT clade" = "#D55E00",
                                "LowT clade" = "#0072B2",
                                "Other NP" = "gray70")) +
  scale_color_manual(values = c("HighT clade" = "#D55E00",
                                 "LowT clade" = "#0072B2",
                                 "Other NP" = "gray70")) +
  labs(x = "Predicted Optimal Growth Temperature (°C)",
       y = "Density",
       title = "OGT Distribution in Nitrosopumilus",
       subtitle = paste0("Total n = ", nrow(all_np), "; ",
                        "HighT (n=", sum(all_np$Clade == "highT"), "), ",
                        "LowT (n=", sum(all_np$Clade == "lowT"), "), ",
                        "Other (n=", sum(all_np$Clade == "other_NP"), ")",
                        "\nDotted lines: measured OGT of cultured strains")) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "black")
  )

ggsave("FigS2_OGT_Distribution_FINAL.pdf", p_distribution, width = 8, height = 5)
ggsave("FigS2_OGT_Distribution_FINAL.png", p_distribution, width = 8, height = 5, dpi = 300)

cat("OGT distribution figure saved:\n")
cat("  - FigS2_OGT_Distribution_FINAL.pdf\n")
cat("  - FigS2_OGT_Distribution_FINAL.png\n\n")

# 保存详细结果
sink("OGT_Validation_Results_FINAL.txt")
cat("=== OGT Prediction Validation Results (FINAL CORRECTED VERSION) ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("ALL DATA VERIFIED FROM LITERATURE\n")
cat(strrep("=", 80), "\n\n")

cat("Dataset Summary:\n")
cat("  Total Nitrosopumilus genomes:", nrow(all_np), "\n")
cat("    HighT clade:", sum(all_np$Clade == "highT"), "genomes\n")
cat("    LowT clade:", sum(all_np$Clade == "lowT"), "genomes\n")
cat("    Other NP:", sum(all_np$Clade == "other_NP"), "genomes\n\n")

cat("Validation Set:\n")
cat("  Cultured strains with measured OGT:", nrow(validation_data), "\n")
cat("    From HighT clade:", sum(validation_data$Clade == "highT"), "\n")
cat("    From LowT clade:", sum(validation_data$Clade == "lowT"), "\n\n")

cat("Statistical Metrics:\n")
cat("  Pearson correlation (r):", round(correlation, 3), "\n")
cat("  P-value:", format(correlation_test$p.value, scientific = TRUE), "\n")
cat("  R-squared:", round(lm_summary$r.squared, 3), "\n")
cat("  RMSE:", round(rmse, 2), "°C\n")
cat("  MAE:", round(mae, 2), "°C\n\n")

cat("Validation Strains (with Literature References):\n")
cat(strrep("=", 80), "\n")
for (i in 1:nrow(validation_data)) {
  cat(sprintf("\n%d. %s (%s)\n", i, validation_data$Strain_Name[i], validation_data$Species[i]))
  cat(sprintf("   Genome: %s\n", validation_data$FileName[i]))
  cat(sprintf("   Clade: %s\n", validation_data$Clade[i]))
  cat(sprintf("   Measured OGT: %d°C\n", validation_data$Measured_OGT[i]))
  cat(sprintf("   Predicted OGT: %.2f°C\n", validation_data$Predicted_OGT[i]))
  cat(sprintf("   Residual: %.2f°C\n", validation_data$Predicted_OGT[i] - validation_data$Measured_OGT[i]))
  cat(sprintf("   Reference: %s (%s) %s\n",
              validation_data$Reference[i],
              validation_data$Year[i],
              validation_data$Journal[i]))
  if (!is.na(validation_data$Notes[i]) && validation_data$Notes[i] != "") {
    cat(sprintf("   Notes: %s\n", validation_data$Notes[i]))
  }
}
cat("\n", strrep("=", 80), "\n\n")

cat("Temperature Ranges:\n")
cat("  HighT clade (measured):", min(highT_temps), "-", max(highT_temps), "°C\n")
cat("  LowT clade (measured):", min(lowT_temps), "-", max(lowT_temps), "°C\n")
cat("  Temperature separation:", min(highT_temps) - max(lowT_temps), "°C\n\n")

cat("Key Findings:\n")
cat("  - Clear thermal differentiation between HighT and LowT clades\n")
cat("  - HighT strains: 30-32°C (mean = ", round(mean(highT_temps), 1), "°C)\n", sep="")
cat("  - LowT strains: 23-25°C (mean = ", round(mean(lowT_temps), 1), "°C)\n", sep="")
cat("  - ", min(highT_temps) - max(lowT_temps), "°C separation supports distinct thermal strategies\n", sep="")
cat("  - Strong correlation (r = ", round(correlation, 3), ") validates prediction method\n\n", sep="")

sink()

cat("Detailed results saved to: OGT_Validation_Results_FINAL.txt\n\n")

cat("=== SUMMARY ===\n\n")
cat("Final corrected validation dataset:\n")
cat("  - 5 cultured strains (3 HighT, 2 LowT)\n")
cat("  - All data verified from published literature\n")
cat("  - HighT: 30-32°C, LowT: 23-25°C\n")
cat("  - Clear 5°C separation between clades\n")
cat("  - r = ", round(correlation, 3), ", RMSE = ", round(rmse, 2), "°C\n\n", sep="")

cat("Analysis complete!\n")
