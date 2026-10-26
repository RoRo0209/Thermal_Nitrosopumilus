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
## Final Analysis: Using ASR-defined clades from chiplot-tree.xlsx
## LowT clade (n=23) and HighT clade (n=11) based on ASR node probabilities
## =========================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(patchwork)
  library(scales)
})

setwd(thermal_path("1.phylogenetic"))

## =========================================================
## 1. Load data with ASR-defined clades
## =========================================================

# Tree
tree_ml <- read.tree("gtdbtk.ar53.NP.treefile")
tree_chrono <- chronos(tree_ml, lambda = 1, quiet = TRUE)
tree <- tree_ml
tree$edge.length <- tree_chrono$edge.length

# Metadata with ASR-defined groups
meta <- read.xlsx("chiplot-tree.xlsx", sheet = 1)
names(meta)[names(meta) == "predOGT.(C)"] <- "OGT"

# The 'group' column contains ASR-defined clades
cat("=== ASR-DEFINED CLADES (from chiplot-tree.xlsx) ===\n")
cat("  LowT clade:", sum(meta$group == "lowT", na.rm = TRUE), "genomes\n")
cat("  HighT clade:", sum(meta$group == "highT", na.rm = TRUE), "genomes\n")
cat("  Not in main clades:", sum(is.na(meta$group)), "genomes\n\n")

# Rename for clarity
meta$Clade <- NA
meta$Clade[meta$group == "lowT"] <- "LowT"
meta$Clade[meta$group == "highT"] <- "HighT"
meta$Clade[is.na(meta$group)] <- "Other"

# Align tree and metadata
common <- intersect(tree$tip.label, meta$FileName)
tree <- drop.tip(tree, setdiff(tree$tip.label, common))
meta <- meta[match(tree$tip.label, meta$FileName), ]

## =========================================================
## 2. OGT statistics for ASR-defined clades
## =========================================================

clade_meta <- meta %>% filter(Clade %in% c("LowT", "HighT"))

cat("=== CLADE OGT SUMMARY ===\n")
ogt_summary <- clade_meta %>%
  group_by(Clade) %>%
  summarise(
    n = n(),
    mean_OGT = mean(OGT, na.rm = TRUE),
    sd_OGT = sd(OGT, na.rm = TRUE),
    min_OGT = min(OGT, na.rm = TRUE),
    max_OGT = max(OGT, na.rm = TRUE)
  )
print(ogt_summary)

t_result <- t.test(OGT ~ Clade, data = clade_meta, var.equal = FALSE)
cat("\nWelch's t-test (LowT clade vs HighT clade):\n")
cat("  t =", round(t_result$statistic, 3), "\n")
cat("  p-value =", format.pval(t_result$p.value, digits = 3), "\n")
cat("  Mean difference =", round(diff(t_result$estimate), 2), "°C\n\n")

## =========================================================
## 3. Phylogenetic signal (Clade members only)
## =========================================================

clade_OGT <- setNames(clade_meta$OGT, clade_meta$FileName)
clade_tree <- keep.tip(tree, names(clade_OGT))

cat("=== PHYLOGENETIC SIGNAL (Clade members only) ===\n")

K_result <- phylosig(clade_tree, clade_OGT, method = "K", test = TRUE, nsim = 1000)
cat("Blomberg's K:\n")
cat("  K =", round(K_result$K, 4), "\n")
cat("  p-value =", format.pval(K_result$P, digits = 3), "\n\n")

lambda_result <- phylosig(clade_tree, clade_OGT, method = "lambda", test = TRUE)
cat("Pagel's lambda:\n")
cat("  λ =", round(lambda_result$lambda, 4), "\n")
cat("  p-value =", format.pval(lambda_result$P, digits = 3), "\n\n")

## =========================================================
## 4. Visualization
## =========================================================

# Color scheme
color_scheme <- c(
  "LowT"  = "#4575B4",
  "HighT" = "#D73027",
  "Other" = "grey80"
)

# Panel A: Circular tree highlighting clades
p_tree <- ggtree(tree, layout = "circular", size = 0.3) %<+% meta +
  geom_tippoint(aes(color = Clade, size = Clade), alpha = 0.8) +
  scale_color_manual(values = color_scheme,
                     name = "Clade",
                     labels = c("LowT" = "LowT clade (ASR)",
                                "HighT" = "HighT clade (ASR)",
                                "Other" = "Other")) +
  scale_size_manual(values = c("LowT" = 3, "HighT" = 3, "Other" = 1),
                    guide = "none") +
  theme(legend.position = c(0.1, 0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10))

# Panel B: OGT distribution
pB <- ggplot(clade_meta, aes(x = OGT, fill = Clade)) +
  geom_histogram(aes(y = after_stat(density)), bins = 12,
                 alpha = 0.7, position = "identity") +
  geom_density(aes(color = Clade), linewidth = 1.2, fill = NA) +
  scale_fill_manual(values = color_scheme) +
  scale_color_manual(values = color_scheme) +
  labs(x = "Optimal Growth Temperature (°C)",
       y = "Density",
       title = "OGT Distribution (ASR-defined Clades)") +
  theme_classic(base_size = 11) +
  theme(legend.position = c(0.8, 0.8),
        legend.background = element_rect(fill = "white", color = "black"),
        plot.title = element_text(face = "bold"))

# Panel C: Phylogenetic signal
signal_data <- data.frame(
  Metric = c("Blomberg's K", "Pagel's λ"),
  Value = c(K_result$K, lambda_result$lambda),
  P_value = c(K_result$P, lambda_result$P),
  Significant = c(K_result$P < 0.05, lambda_result$P < 0.05)
)

pC <- ggplot(signal_data, aes(x = Metric, y = Value, fill = Significant)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0("p = ", format.pval(P_value, digits = 2))),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#2166AC")) +
  labs(y = "Statistic Value", x = NULL,
       title = "Phylogenetic Signal") +
  ylim(0, max(signal_data$Value) * 1.2) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

# Panel D: Boxplot
pD <- ggplot(clade_meta, aes(x = Clade, y = OGT, fill = Clade)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2.5) +
  scale_fill_manual(values = color_scheme) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 4,
               fill = "white", color = "black", stroke = 1.5) +
  labs(y = "OGT (°C)", x = NULL,
       title = "OGT Comparison") +
  annotate("text", x = 1.5, y = max(clade_meta$OGT, na.rm = TRUE) + 1,
           label = paste0("ΔT = ", round(abs(diff(t_result$estimate)), 1), "°C\n",
                         "p < ", format.pval(t_result$p.value, digits = 1)),
           size = 4, fontface = "bold") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

## =========================================================
## 5. Combine and save
## =========================================================

combined <- (p_tree | (pB / pC)) +
  plot_annotation(tag_levels = "A",
                  title = "Phylogenetic Framework of AOA Ecotypes (ASR-defined Clades)",
                  theme = theme(plot.title = element_text(face = "bold", size = 14),
                                plot.tag = element_text(face = "bold", size = 13)))

ggsave("Figure1_ASR_clades_final.pdf", combined,
       width = 16, height = 10, dpi = 300, bg = "white")

ggsave("Figure1_ASR_clades_final.png", combined,
       width = 16, height = 10, dpi = 300, bg = "white")

# Save tree only
ggsave("Phylo_tree_ASR_clades.pdf", p_tree,
       width = 10, height = 10, bg = "white")

# Save clade assignments
write.csv(meta %>% select(FileName, OGT, group, Clade),
          "ASR_clade_assignments_final.csv", row.names = FALSE)

# Save clade member lists for downstream analysis
lowT_members <- meta$FileName[meta$Clade == "LowT"]
highT_members <- meta$FileName[meta$Clade == "HighT"]

writeLines(lowT_members, "LowT_clade_final.txt")
writeLines(highT_members, "HighT_clade_final.txt")

cat("\n=== OUTPUT ===\n")
cat("Figures saved:\n")
cat("  - Figure1_ASR_clades_final.pdf/png (main figure)\n")
cat("  - Phylo_tree_ASR_clades.pdf (tree only)\n")
cat("  - ASR_clade_assignments_final.csv (all genomes)\n")
cat("  - LowT_clade_final.txt (", length(lowT_members), " genomes)\n")
cat("  - HighT_clade_final.txt (", length(highT_members), " genomes)\n\n")

cat("=== FINAL SUMMARY ===\n")
cat("Supplied clade assignments (not inferred by this script):\n")
cat("  LowT clade:", length(lowT_members), "genomes, mean OGT =",
    round(mean(clade_meta$OGT[clade_meta$Clade == "LowT"], na.rm = TRUE), 2), "°C\n")
cat("  HighT clade:", length(highT_members), "genomes, mean OGT =",
    round(mean(clade_meta$OGT[clade_meta$Clade == "HighT"], na.rm = TRUE), 2), "°C\n")
cat("  ΔT =", round(abs(diff(t_result$estimate)), 2), "°C (p =",
    format.pval(t_result$p.value, digits = 3), ")\n")
cat("  Phylogenetic signal: K =", round(K_result$K, 3),
    "(p =", format.pval(K_result$P, digits = 3), "),",
    "λ =", round(lambda_result$lambda, 3),
    "(p =", format.pval(lambda_result$P, digits = 3), ")\n\n")

