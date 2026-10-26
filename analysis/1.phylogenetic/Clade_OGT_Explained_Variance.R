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

## =========================================================
## Clade membership explains OGT variation
## Linear model: OGT ~ Clade, reporting R², η², Cohen's d
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
})

## Auto-detect Windows vs WSL path
if (dir.exists(thermal_path("1.phylogenetic"))) {
  setwd(thermal_path("1.phylogenetic"))
} else {
  setwd(thermal_path("1.phylogenetic"))
}

## 1. Load data ----
dat <- read.csv("ASR_clade_assignments_final.csv", stringsAsFactors = FALSE)
clade_dat <- dat %>% filter(Clade %in% c("LowT", "HighT"))
clade_dat$Clade <- factor(clade_dat$Clade, levels = c("LowT", "HighT"))

cat("=== Sample sizes ===\n")
cat("  LowT:", sum(clade_dat$Clade == "LowT"), "\n")
cat("  HighT:", sum(clade_dat$Clade == "HighT"), "\n\n")

## 2. Linear model: OGT ~ Clade ----
fit <- lm(OGT ~ Clade, data = clade_dat)
cat("=== Linear Model: OGT ~ Clade ===\n")
print(summary(fit))

r2 <- summary(fit)$r.squared
adj_r2 <- summary(fit)$adj.r.squared
cat("\n  R²     =", round(r2, 4), "\n")
cat("  Adj R² =", round(adj_r2, 4), "\n")

## 3. ANOVA → η² (eta-squared) ----
aov_fit <- aov(OGT ~ Clade, data = clade_dat)
aov_tab <- summary(aov_fit)[[1]]
SS_clade <- aov_tab["Clade", "Sum Sq"]
SS_total <- sum(aov_tab[, "Sum Sq"])
eta_sq <- SS_clade / SS_total

cat("\n=== ANOVA ===\n")
print(aov_tab)
cat("\n  η² (eta-squared) =", round(eta_sq, 4), "\n")
cat("  (Same as R² for single-predictor model)\n")

## 4. Cohen's d ----
lowT_ogt <- clade_dat$OGT[clade_dat$Clade == "LowT"]
highT_ogt <- clade_dat$OGT[clade_dat$Clade == "HighT"]

n1 <- length(lowT_ogt)
n2 <- length(highT_ogt)
pooled_sd <- sqrt(((n1 - 1) * sd(lowT_ogt)^2 + (n2 - 1) * sd(highT_ogt)^2) / (n1 + n2 - 2))
cohens_d <- abs(mean(highT_ogt) - mean(lowT_ogt)) / pooled_sd

cat("\n=== Effect Size ===\n")
cat("  Cohen's d =", round(cohens_d, 3), "\n")
cat("  Pooled SD =", round(pooled_sd, 3), "°C\n")
cat("  Interpretation: d > 0.8 = large effect\n\n")

## 5. Summary for manuscript ----
cat("=== FOR MANUSCRIPT ===\n")
cat(sprintf(
  "Clade identity explained %.1f%% of OGT variation among the 34 genomes\n",
  r2 * 100
))
cat(sprintf(
  "(linear model: R² = %.3f, F(1,%d) = %.1f, p < 2×10⁻¹⁶; η² = %.3f; Cohen's d = %.2f).\n",
  r2, n1 + n2 - 2,
  aov_tab["Clade", "F value"],
  eta_sq, cohens_d
))

cat("\nDONE\n")
