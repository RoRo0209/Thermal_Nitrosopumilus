# Release path adapter: data and outputs live in a separate user workspace.
.release_root <- Sys.getenv("THERMAL_PROJECT_ROOT")
if (!nzchar(.release_root)) stop("Set THERMAL_PROJECT_ROOT to the analysis workspace")
.release_root <- normalizePath(.release_root, winslash = "/", mustWork = TRUE)
.release_input <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) stop("Set required input environment variable: ", name)
  normalizePath(value, winslash = "/", mustWork = TRUE)
}
setwd(file.path(.release_root, "instrain/out"))

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(rstatix)
})

setwd(file.path(.release_root, "instrain", "out"))

df0 <- read_excel("instrain_coverage_dayu1.xlsx", sheet = 1)

df <- df0 %>%
  transmute(
    group = factor(clade, levels = c("lowT", "highT")),
    nucleotide_diversity = as.numeric(nucl_diversity_rarefied),
    SNV_per_kbp         = as.numeric(SNV_count) / as.numeric(length) * 1000,
    pN_pS               = as.numeric(pN_pS),
    D_prime             = as.numeric(d_prime_mean)
  ) %>%
  filter(!is.na(group))

df_long <- df %>%
  pivot_longer(-group, names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))

metric_labs <- c(
  nucleotide_diversity = "Nucleotide diversity (π)",
  SNV_per_kbp          = "SNVs per kbp",
  pN_pS                = "pN/pS",
  D_prime              = "D'"
)

group_cols <- c(lowT = "#8EC6E8", highT = "#F4A261")

# Wilcoxon test per metric
stat_res <- df_long %>%
  group_by(metric) %>%
  wilcox_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(p_label = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01  ~ "**",
    p.adj <= 0.05  ~ "*",
    TRUE ~ "ns"
  )) %>%
  left_join(
    df_long %>% group_by(metric) %>%
      summarise(y_pos = max(value, na.rm = TRUE) * 1.08, .groups = "drop"),
    by = "metric"
  )

print(stat_res %>% select(metric, p, p.adj, p_label))

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor    = element_blank(),
    legend.position     = "top"
  )

make_panel <- function(m) {
  dat <- df_long %>% filter(metric == m)
  st  <- stat_res %>% filter(metric == m)
  ggplot(dat, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.5, show.legend = FALSE) +
    geom_text(data = st, aes(x = 1.5, y = y_pos, label = p_label),
              inherit.aes = FALSE, size = 4) +
    scale_fill_manual(values = group_cols) +
    labs(x = NULL, y = metric_labs[[m]], fill = "Group") +
    base_theme
}

fig <- (make_panel("nucleotide_diversity") | make_panel("SNV_per_kbp")) /
       (make_panel("pN_pS")               | make_panel("D_prime")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave("instrain_boxplot.png", fig, width = 8, height = 7, dpi = 300)
ggsave("instrain_boxplot.pdf", fig, width = 8, height = 7)
message("Done.")
