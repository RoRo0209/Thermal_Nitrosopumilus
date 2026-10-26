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
# 黑潮–亲潮过渡带 NP 生态型空间梯度验证（v3）
# 修复: (1) PERMANOVA by="terms"; (2) 新增环境变量相关性分析(图04)
# ============================================================================

suppressPackageStartupMessages({
  library(readxl); library(tibble); library(dplyr); library(tidyr)
  library(stringr); library(forcats); library(purrr); library(vegan)
  library(ggplot2); library(ggrepel); library(cowplot)
  library(patchwork); library(scales)
})

# --- 路径设置 ---
script_dir <- file.path(.release_root, "kuroshio")
setwd(script_dir)
dir.create("output_figures", showWarnings = FALSE)
dir.create("output_tables",  showWarnings = FALSE)

# --- 颜色方案 ---
COL <- list(
  highT      = "#D4541A",
  lowT       = "#1E6FA4",
  transition = "#7A9E7E",
  kuroshio   = "#E8913A",
  oyashio    = "#3A7EC2",
  mixed      = "#8CAE8C"
)

# --- 用户选项 ---
wm_mode                <- "TS_threshold"   # "TS_threshold" or "TS_kmeans"
include_KEmixed_in_warm <- FALSE

# --- T-S 阈值（文献参考值，可调整）---
thr <- list(
  S_kw = 34.4, T_kw = 20,
  T_kem_lo = 15, T_kem_hi = 20,
  T_oy = 10, S_oy = 34.0, T_oy2 = 15
)

# ============================================================================
# 0) 读取数据
# ============================================================================

# 0.1 环境因子
env_raw <- readxl::read_excel("environmental_factor.xlsx", sheet = 1, skip = 1) %>%
  as.data.frame()

colnames(env_raw) <- c(
  "Station", "Salinity", "Temperature", "Pressure",
  "DO", "Chla", "Ammonium", "Nitrate", "Nitrite",
  "Phosphate", "Silicate", "Density", "CDOM"
)

env <- env_raw %>%
  dplyr::filter(grepl("^20\\d\\d_D", Station)) %>%
  dplyr::mutate(dplyr::across(-Station, ~ suppressWarnings(as.numeric(.)))) %>%
  dplyr::mutate(
    Year    = as.integer(substr(Station, 1, 4)),
    Section = stringr::str_extract(Station, "D\\d+"),
    D_num   = suppressWarnings(as.integer(stringr::str_extract(Section, "\\d+"))),
    Depth   = stringr::str_extract(Station, "\\d+m$|DCM$"),
    Depth_m = dplyr::case_when(
      Depth == "DCM" ~ 35.3,
      is.na(Depth)   ~ NA_real_,
      TRUE           ~ suppressWarnings(as.numeric(stringr::str_remove(Depth, "m")))
    )
  )

cat("env rows:", nrow(env), "\n")

# 0.2 丰度表
abund_raw <- readxl::read_excel("kuroshio_np_temp_abudance.xlsx", sheet = 1) %>%
  as.data.frame()

genome_col <- names(abund_raw)[1]
group_col  <- names(abund_raw)[stringr::str_detect(tolower(names(abund_raw)), "^group$|ecotype|type|cluster")][1]
if (is.na(group_col)) stop("Cannot find group/ecotype/type column in abundance table.")

abund_long <- abund_raw %>%
  dplyr::rename(Genome = !!genome_col, Group = !!group_col) %>%
  tidyr::pivot_longer(cols = -c(Genome, Group), names_to = "Station", values_to = "Abundance") %>%
  dplyr::mutate(
    Abundance = suppressWarnings(as.numeric(as.character(Abundance))),
    Abundance = tidyr::replace_na(Abundance, 0),
    Group     = tolower(trimws(as.character(Group)))
  )

abund_wide <- abund_long %>%
  dplyr::select(Station, Genome, Abundance) %>%
  tidyr::pivot_wider(names_from = Genome, values_from = Abundance, values_fill = 0) %>%
  tibble::column_to_rownames("Station")

cat("abundance matrix:", nrow(abund_wide), "stations x", ncol(abund_wide), "genomes\n")

# 0.3 汇总 HighT/LowT 并合并
abund_summary <- abund_long %>%
  dplyr::group_by(Station, Group) %>%
  dplyr::summarise(Total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Group, values_from = Total, values_fill = 0) %>%
  dplyr::rename_with(~ dplyr::case_when(
    tolower(.) == "hight" ~ "HighT",
    tolower(.) == "lowt"  ~ "LowT",
    TRUE ~ .
  )) %>%
  dplyr::mutate(Sum = HighT + LowT)

common_sta <- intersect(env$Station, abund_summary$Station)

meta <- env %>%
  dplyr::filter(Station %in% common_sta) %>%
  dplyr::inner_join(abund_summary %>% dplyr::filter(Station %in% common_sta), by = "Station") %>%
  dplyr::filter(is.finite(Temperature) & is.finite(Salinity)) %>%
  dplyr::filter(Sum > 0) %>%
  dplyr::mutate(
    HighT_rel = HighT / Sum,
    LowT_rel  = LowT  / Sum,
    Dominance = dplyr::case_when(
      HighT_rel >= 0.6 ~ "HighT-dominated",
      LowT_rel  >= 0.6 ~ "LowT-dominated",
      TRUE             ~ "Mixed"
    )
  )

write.csv(meta, "output_tables/meta_with_env_and_sum.csv", row.names = FALSE)
cat("meta stations (Sum>0):", nrow(meta), "\n")

# ============================================================================
# 1) 水团分类
# ============================================================================
if (wm_mode == "TS_threshold") {
  meta <- meta %>%
    dplyr::mutate(
      WaterMass = dplyr::case_when(
        Temperature >= thr$T_kw  & Salinity >= thr$S_kw                          ~ "Kuroshio",
        Temperature >= thr$T_kem_lo & Temperature < thr$T_kem_hi & Salinity >= thr$S_kw ~ "KE_Mixed",
        Temperature <  thr$T_oy                                                   ~ "Oyashio",
        Temperature <  thr$T_oy2 & Salinity < thr$S_oy                           ~ "Oyashio",
        TRUE ~ "Transition"
      ),
      WaterMass = factor(WaterMass, levels = c("Kuroshio","KE_Mixed","Transition","Oyashio"))
    )
} else if (wm_mode == "TS_kmeans") {
  ts_sc <- scale(meta %>% dplyr::select(Temperature, Salinity) %>% as.matrix())
  set.seed(1)
  km <- stats::kmeans(ts_sc, centers = 3, nstart = 50)
  meta$WM_km <- paste0("C", km$cluster)
  cl_means <- meta %>%
    dplyr::group_by(WM_km) %>%
    dplyr::summarise(Tm = mean(Temperature), Sm = mean(Salinity), .groups = "drop")
  warm_id <- cl_means$WM_km[which.max(cl_means$Tm + cl_means$Sm)]
  cold_id <- cl_means$WM_km[which.min(cl_means$Tm + cl_means$Sm)]
  meta <- meta %>%
    dplyr::mutate(
      WaterMass = dplyr::case_when(
        WM_km == warm_id ~ "Kuroshio",
        WM_km == cold_id ~ "Oyashio",
        TRUE ~ "Transition"
      ),
      WaterMass = factor(WaterMass, levels = c("Kuroshio","Transition","Oyashio"))
    )
} else {
  stop("wm_mode must be 'TS_threshold' or 'TS_kmeans'")
}

cat("\nWater mass counts (mode =", wm_mode, "):\n")
print(table(meta$WaterMass))
write.csv(meta %>% dplyr::select(Station, Temperature, Salinity, WaterMass),
          "output_tables/watermass_labels.csv", row.names = FALSE)

# --- 图01: T-S 图 ---
p_ts <- ggplot2::ggplot(meta, ggplot2::aes(x = Salinity, y = Temperature)) +
  ggplot2::geom_point(ggplot2::aes(colour = WaterMass, size = HighT_rel + 0.02),
                      alpha = 0.85, stroke = 0.5) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = stringr::str_remove(Station, "^20\\d\\d_")),
    size = 2.4, max.overlaps = 25, colour = "grey35"
  ) +
  ggplot2::scale_colour_manual(
    values = c(Kuroshio = COL$kuroshio, KE_Mixed = "#BF8A3A",
               Transition = COL$transition, Oyashio = COL$oyashio),
    name = "Water mass"
  ) +
  ggplot2::scale_size_continuous(range = c(2, 9), name = "HighT rel. abund.") +
  ggplot2::labs(
    title = "T-S diagram (water mass classification)",
    subtitle = paste0("Mode = ", wm_mode),
    x = "Salinity (PSU)", y = "Temperature (degC)"
  ) +
  cowplot::theme_cowplot(12)

if (wm_mode == "TS_threshold") {
  p_ts <- p_ts +
    ggplot2::geom_hline(yintercept = thr$T_kw, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = thr$S_kw, linetype = "dashed", colour = "grey55", linewidth = 0.4)
}

ggplot2::ggsave("output_figures/01_TS_diagram.pdf", p_ts, width = 9, height = 6)
ggplot2::ggsave("output_figures/01_TS_diagram.png", p_ts, width = 9, height = 6, dpi = 300)

# ============================================================================
# 2) 暖水团 vs 冷水团比较
# ============================================================================
warm_levels <- if (include_KEmixed_in_warm) c("Kuroshio","KE_Mixed") else c("Kuroshio")

meta_kw <- meta %>%
  dplyr::filter(WaterMass %in% c(warm_levels, "Oyashio")) %>%
  dplyr::mutate(
    WM2 = ifelse(WaterMass %in% warm_levels, "Warm (Kuroshio)", "Cold (Oyashio)"),
    WM2 = factor(WM2, levels = c("Warm (Kuroshio)", "Cold (Oyashio)"))
  )

wt_hight <- stats::wilcox.test(HighT_rel ~ WM2, data = meta_kw, exact = FALSE)
wt_lowt  <- stats::wilcox.test(LowT_rel  ~ WM2, data = meta_kw, exact = FALSE)

pval_fmt <- function(p) {
  if (p < 0.001) "p < 0.001 (***)"
  else if (p < 0.01)  sprintf("p = %.3f (**)",  p)
  else if (p < 0.05)  sprintf("p = %.3f (*)",   p)
  else                sprintf("p = %.3f (ns)",   p)
}

make_box <- function(df, yvar, ytitle, wt) {
  ymax <- max(df[[yvar]], na.rm = TRUE)
  ggplot2::ggplot(df, ggplot2::aes(x = WM2, y = .data[[yvar]], fill = WM2)) +
    ggplot2::geom_boxplot(width = 0.48, alpha = 0.75, outlier.shape = 21, outlier.size = 2.5) +
    ggplot2::geom_jitter(width = 0.14, size = 3, alpha = 0.7, shape = 21) +
    ggplot2::scale_fill_manual(
      values = c("Warm (Kuroshio)" = COL$kuroshio, "Cold (Oyashio)" = COL$oyashio),
      guide = "none"
    ) +
    ggplot2::annotate("segment", x = 1, xend = 2, y = ymax * 1.06, yend = ymax * 1.06, linewidth = 0.5) +
    ggplot2::annotate("text", x = 1.5, y = ymax * 1.13, label = pval_fmt(wt$p.value),
                      size = 3.8, fontface = "bold") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(),
                                expand = ggplot2::expansion(mult = c(0.02, 0.20))) +
    ggplot2::labs(x = NULL, y = ytitle) +
    cowplot::theme_cowplot(12)
}

p_box_h <- make_box(meta_kw, "HighT_rel", "HighT relative abundance", wt_hight)
p_box_l <- make_box(meta_kw, "LowT_rel",  "LowT relative abundance",  wt_lowt)

# 站点组成柱状图（按温度排序）
bar_df <- meta %>%
  dplyr::arrange(desc(Temperature)) %>%
  dplyr::mutate(Sta_ord = forcats::fct_inorder(stringr::str_remove(Station, "^20\\d\\d_"))) %>%
  dplyr::select(Sta_ord, Temperature, WaterMass, HighT_rel, LowT_rel) %>%
  tidyr::pivot_longer(cols = c(HighT_rel, LowT_rel), names_to = "Ecotype", values_to = "Rel") %>%
  dplyr::mutate(Ecotype = dplyr::recode(Ecotype, HighT_rel = "HighT", LowT_rel = "LowT"))

temp_labels <- meta %>%
  dplyr::arrange(desc(Temperature)) %>%
  dplyr::mutate(Sta_ord = forcats::fct_inorder(stringr::str_remove(Station, "^20\\d\\d_"))) %>%
  dplyr::select(Sta_ord, Temperature)

p_bar <- ggplot2::ggplot(bar_df, ggplot2::aes(x = Sta_ord, y = Rel, fill = Ecotype)) +
  ggplot2::geom_col(width = 0.78, colour = "white", linewidth = 0.2) +
  ggplot2::geom_text(data = temp_labels,
                     ggplot2::aes(x = Sta_ord, y = 1.04, label = sprintf("%.1f", Temperature)),
                     size = 2.6, inherit.aes = FALSE) +
  ggplot2::scale_fill_manual(values = c(HighT = COL$highT, LowT = COL$lowT)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(),
                              expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::labs(x = "Station (decreasing temperature ->)",
                y = "Relative abundance", fill = "NP ecotype",
                title = "Ecotype composition along transect") +
  cowplot::theme_cowplot(11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))

p_watermass <- (p_box_h | p_box_l) / p_bar +
  patchwork::plot_annotation(
    title = "A. Water mass comparison: warm vs cold",
    subtitle = paste0("wm_mode=", wm_mode, "; KE_Mixed in warm=", include_KEmixed_in_warm),
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
  )

ggplot2::ggsave("output_figures/02_watermass_comparison.pdf", p_watermass, width = 12, height = 10)
ggplot2::ggsave("output_figures/02_watermass_comparison.png", p_watermass, width = 12, height = 10, dpi = 300)

# ============================================================================
# 3) 生态型替代曲线（沿温度梯度）
# ============================================================================
meta_sc <- meta %>% dplyr::filter(is.finite(Temperature))

cor_h <- stats::cor.test(meta_sc$Temperature, meta_sc$HighT_rel, method = "spearman")
cor_l <- stats::cor.test(meta_sc$Temperature, meta_sc$LowT_rel,  method = "spearman")

long_sc <- meta_sc %>%
  tidyr::pivot_longer(cols = c(HighT_rel, LowT_rel), names_to = "Ecotype", values_to = "Rel") %>%
  dplyr::mutate(Ecotype = dplyr::recode(Ecotype, HighT_rel = "HighT", LowT_rel = "LowT"))

annot_corr <- data.frame(
  Ecotype = c("HighT", "LowT"),
  xpos    = rep(max(meta_sc$Temperature, na.rm = TRUE), 2),
  ypos    = c(0.93, 0.82),
  label   = c(
    sprintf("rho = %+.2f, p = %.3f", cor_h$estimate, cor_h$p.value),
    sprintf("rho = %+.2f, p = %.3f", cor_l$estimate, cor_l$p.value)
  )
)

p_curve <- ggplot2::ggplot(long_sc, ggplot2::aes(x = Temperature, y = Rel, colour = Ecotype)) +
  ggplot2::geom_smooth(method = "loess", span = 0.85, se = TRUE, alpha = 0.15, linewidth = 1.3) +
  ggplot2::geom_point(ggplot2::aes(shape = WaterMass), size = 3.2, alpha = 0.85, stroke = 0.5) +
  ggplot2::geom_text(data = annot_corr,
                     ggplot2::aes(x = xpos, y = ypos, label = label, colour = Ecotype),
                     size = 3.6, fontface = "bold", inherit.aes = FALSE, hjust = 1) +
  ggplot2::scale_colour_manual(values = c(HighT = COL$highT, LowT = COL$lowT)) +
  ggplot2::scale_shape_manual(
    values = c(Kuroshio = 16, KE_Mixed = 17, Transition = 15, Oyashio = 18),
    name = "Water mass"
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.05)) +
  ggplot2::labs(x = "Temperature (degC)", y = "Relative abundance",
                colour = "NP ecotype",
                title = "B. Ecotype replacement along temperature gradient",
                subtitle = "LOESS smoothing (span=0.85)") +
  cowplot::theme_cowplot(12)

ggplot2::ggsave("output_figures/03_replacement_curve.pdf", p_curve, width = 10, height = 6)
ggplot2::ggsave("output_figures/03_replacement_curve.png", p_curve, width = 10, height = 6, dpi = 300)

# ============================================================================
# 4) 【新增】环境变量与 HighT_rel 的相关性分析
# ============================================================================
env_vars <- c("Temperature","Salinity","DO","Chla","Nitrate","Nitrite","Phosphate","Silicate")

cor_env <- purrr::map_dfr(env_vars, function(v) {
  x <- meta[[v]]
  y <- meta$HighT_rel
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(NULL)
  ct <- stats::cor.test(x[ok], y[ok], method = "spearman")
  data.frame(
    Variable = v,
    rho      = as.numeric(ct$estimate),
    p        = ct$p.value,
    n        = sum(ok)
  )
})

cor_env <- cor_env %>%
  dplyr::mutate(
    sig   = dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"),
    label = sprintf("rho=%.2f%s", rho, sig),
    Variable = factor(Variable, levels = env_vars[order(cor_env$rho[match(env_vars, cor_env$Variable)])])
  )

write.csv(cor_env, "output_tables/env_correlation_HighT.csv", row.names = FALSE)

p_cor_env <- ggplot2::ggplot(cor_env, ggplot2::aes(x = rho, y = Variable, fill = rho)) +
  ggplot2::geom_col(width = 0.65, colour = "white") +
  ggplot2::geom_text(ggplot2::aes(
    label = label,
    hjust = ifelse(rho >= 0, -0.1, 1.1)
  ), size = 3.5) +
  ggplot2::scale_fill_gradient2(
    low = COL$oyashio, mid = "grey90", high = COL$kuroshio,
    midpoint = 0, limits = c(-1, 1), name = "Spearman rho"
  ) +
  ggplot2::scale_x_continuous(limits = c(-1.15, 1.15), breaks = seq(-1, 1, 0.5)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
  ggplot2::labs(
    x = "Spearman rho (vs HighT relative abundance)",
    y = NULL,
    title = "C. Environmental correlates of HighT ecotype dominance",
    subtitle = "Spearman correlation; * p<0.05, ** p<0.01, *** p<0.001"
  ) +
  cowplot::theme_cowplot(12) +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave("output_figures/04_env_correlation.pdf", p_cor_env, width = 8, height = 5)
ggplot2::ggsave("output_figures/04_env_correlation.png", p_cor_env, width = 8, height = 5, dpi = 300)

# ============================================================================
# 5) Mantel + Partial Mantel
# ============================================================================
common_sta2 <- intersect(rownames(abund_wide), meta$Station)
abund_m <- abund_wide[common_sta2, , drop = FALSE]

env_m <- meta %>%
  dplyr::filter(Station %in% common_sta2) %>%
  dplyr::arrange(match(Station, common_sta2)) %>%
  dplyr::select(Station, Temperature, Salinity, DO, Chla, Nitrate, Phosphate) %>%
  tibble::column_to_rownames("Station")

abund_rel_m <- vegan::decostand(abund_m, method = "total")
dist_bc     <- vegan::vegdist(abund_rel_m, method = "bray")

tvec <- env_m$Temperature; svec <- env_m$Salinity
tvec[!is.finite(tvec)] <- median(tvec, na.rm = TRUE)
svec[!is.finite(svec)] <- median(svec, na.rm = TRUE)

dist_temp <- as.dist(abs(outer(tvec, tvec, "-")))
dist_sal  <- as.dist(abs(outer(svec, svec, "-")))

env_num <- env_m %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(!is.finite(.), median(., na.rm = TRUE), .)))
dist_env <- stats::dist(scale(env_num), method = "euclidean")

set.seed(42)
nperm <- 9999
m_temp <- vegan::mantel(dist_bc, dist_temp, method = "spearman", permutations = nperm, na.rm = TRUE)
m_sal  <- vegan::mantel(dist_bc, dist_sal,  method = "spearman", permutations = nperm, na.rm = TRUE)
m_env  <- vegan::mantel(dist_bc, dist_env,  method = "spearman", permutations = nperm, na.rm = TRUE)
pm_t   <- vegan::mantel.partial(dist_bc, dist_temp, dist_sal, method = "spearman", permutations = nperm, na.rm = TRUE)

mantel_tbl <- data.frame(
  Test = c("Mantel: BC ~ abs(dTemp)", "Mantel: BC ~ abs(dSal)",
           "Mantel: BC ~ Env(all)", "Partial: BC ~ abs(dTemp) | abs(dSal)"),
  r = c(m_temp$statistic, m_sal$statistic, m_env$statistic, pm_t$statistic),
  p = c(m_temp$signif,    m_sal$signif,    m_env$signif,    pm_t$signif)
)
write.csv(mantel_tbl, "output_tables/mantel_results.csv", row.names = FALSE)

dist_df <- data.frame(BC = as.vector(dist_bc), dTemp = as.vector(dist_temp), dSal = as.vector(dist_sal))
write.csv(dist_df, "output_tables/pairwise_distances.csv", row.names = FALSE)

p_m1 <- ggplot2::ggplot(dist_df, ggplot2::aes(x = dTemp, y = BC)) +
  ggplot2::geom_point(alpha = 0.30, size = 1.5, colour = COL$highT) +
  ggplot2::geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  ggplot2::annotate("text",
                    x = max(dist_df$dTemp) * 0.62, y = max(dist_df$BC) * 0.12,
                    label = sprintf("r = %.3f\np = %.4f", m_temp$statistic, m_temp$signif),
                    size = 4, hjust = 0) +
  ggplot2::labs(x = "abs(dTemperature) (degC)", y = "Bray-Curtis dissimilarity",
                title = "Mantel: community ~ temperature difference") +
  cowplot::theme_cowplot(12)

p_m2 <- ggplot2::ggplot(dist_df, ggplot2::aes(x = dSal, y = BC)) +
  ggplot2::geom_point(alpha = 0.30, size = 1.5, colour = COL$lowT) +
  ggplot2::geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  ggplot2::annotate("text",
                    x = max(dist_df$dSal) * 0.62, y = max(dist_df$BC) * 0.12,
                    label = sprintf("r = %.3f\np = %.4f", m_sal$statistic, m_sal$signif),
                    size = 4, hjust = 0) +
  ggplot2::labs(x = "abs(dSalinity) (PSU)", y = "Bray-Curtis dissimilarity",
                title = "Mantel: community ~ salinity difference") +
  cowplot::theme_cowplot(12)

p_mantel <- (p_m1 | p_m2) +
  patchwork::plot_annotation(
    title = "D. Mantel tests (pairwise distances)",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
  )

ggplot2::ggsave("output_figures/05_mantel_scatter.pdf", p_mantel, width = 12, height = 5)
ggplot2::ggsave("output_figures/05_mantel_scatter.png", p_mantel, width = 12, height = 5, dpi = 300)

# ============================================================================
# 6) NMDS + PERMANOVA（修复：by="terms" 获取各变量独立 R2）
# ============================================================================
set.seed(42)
nmds <- vegan::metaMDS(abund_rel_m, distance = "bray", k = 2, trymax = 100, trace = FALSE)

nmds_df <- as.data.frame(vegan::scores(nmds, display = "sites")) %>%
  tibble::rownames_to_column("Station") %>%
  dplyr::left_join(meta %>% dplyr::select(Station, Temperature, WaterMass, HighT_rel, LowT_rel),
                   by = "Station")

perm_env <- env_m[, c("Temperature","Salinity"), drop = FALSE]
ok <- stats::complete.cases(perm_env)

# 关键修复：by = "terms" 才能获得各变量独立的 R2
perm_res <- vegan::adonis2(
  abund_rel_m[ok, , drop = FALSE] ~ Temperature + Salinity,
  data         = perm_env[ok, , drop = FALSE],
  method       = "bray",
  permutations = 9999,
  by           = "terms"   # FIX: was missing, caused empty r2_temp/r2_sal
)

perm_tbl <- as.data.frame(perm_res) %>% tibble::rownames_to_column("Term")
write.csv(perm_tbl, "output_tables/permanova_results.csv", row.names = FALSE)

r2_temp <- perm_tbl$R2[perm_tbl$Term == "Temperature"]
r2_sal  <- perm_tbl$R2[perm_tbl$Term == "Salinity"]

p_nmds <- ggplot2::ggplot(nmds_df, ggplot2::aes(x = NMDS1, y = NMDS2)) +
  ggplot2::geom_point(ggplot2::aes(colour = Temperature, size = HighT_rel + 0.02, shape = WaterMass),
                      alpha = 0.85, stroke = 0.4) +
  ggplot2::scale_colour_gradient2(
    low = COL$oyashio, mid = COL$transition, high = COL$kuroshio,
    midpoint = 17, name = "Temp (degC)"
  ) +
  ggplot2::scale_size_continuous(range = c(3, 9), name = "HighT rel.") +
  ggplot2::scale_shape_manual(values = c(Kuroshio = 16, KE_Mixed = 17, Transition = 15, Oyashio = 18)) +
  ggrepel::geom_text_repel(ggplot2::aes(label = stringr::str_remove(Station, "^20\\d\\d_")),
                           size = 2.5, max.overlaps = 20, colour = "grey30") +
  ggplot2::annotate("text",
                    x = min(nmds_df$NMDS1, na.rm = TRUE), y = max(nmds_df$NMDS2, na.rm = TRUE),
                    label = paste0("Stress = ", round(nmds$stress, 3)),
                    size = 3.5, hjust = 0, fontface = "italic") +
  ggplot2::labs(title = "E. NMDS (Bray-Curtis)",
                subtitle = "Colour = temperature; size = HighT dominance") +
  cowplot::theme_cowplot(12)

ggplot2::ggsave("output_figures/06_NMDS.pdf", p_nmds, width = 9, height = 6)
ggplot2::ggsave("output_figures/06_NMDS.png", p_nmds, width = 9, height = 6, dpi = 300)

# ============================================================================
# 7) 汇总面板
# ============================================================================
p_summary <- (p_ts | p_curve) / (p_m1 | p_nmds) +
  patchwork::plot_annotation(
    title    = "Kuroshio-Oyashio transect: spatial gradient verification",
    subtitle = paste0(
      nrow(meta), " stations (Sum>0); wm_mode=", wm_mode,
      "; PERMANOVA R2(Temp)=", sprintf("%.3f", r2_temp),
      ", R2(Sal)=",  sprintf("%.3f", r2_sal)
    ),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "grey40")
    )
  ) +
  patchwork::plot_layout(guides = "collect")

ggplot2::ggsave("output_figures/00_summary_4panel.pdf", p_summary, width = 16, height = 12)
ggplot2::ggsave("output_figures/00_summary_4panel.png", p_summary, width = 16, height = 12, dpi = 300)

# ============================================================================
# 8) 关键结果汇总
# ============================================================================
cat("\n==================== KEY SUMMARY ====================\n")
cat("wm_mode:", wm_mode, "\n")
cat("WaterMass counts:\n"); print(table(meta$WaterMass))
cat(sprintf("Wilcoxon HighT (warm vs cold): p = %.4f\n", wt_hight$p.value))
cat(sprintf("Spearman HighT ~ Temp: rho = %+.3f, p = %.4f\n", cor_h$estimate, cor_h$p.value))
cat(sprintf("Mantel BC ~ abs(dTemp): r = %.4f, p = %.4f\n", m_temp$statistic, m_temp$signif))
cat(sprintf("Partial Mantel (dTemp|dSal): r = %.4f, p = %.4f\n", pm_t$statistic, pm_t$signif))
cat(sprintf("PERMANOVA R2(Temp) = %.3f; R2(Sal) = %.3f\n", r2_temp, r2_sal))
cat("Outputs: output_figures/ and output_tables/\n")
cat("=====================================================\n")
