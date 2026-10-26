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
setwd(file.path(.release_root, "60weeks"))

# ============================================================
# 60 weeks metagenome points ONLY (E/W/Z split)
# ------------------------------------------------------------
# Goal:
#   - Automatically split samples by station (E/W/Z)
#   - KEEP ONLY metagenome-measured points:
#       i.e., weeks where (HighT + LowT) > 0
#   - For each station: make plots + run GAM/CCF/segmented/logistic
#     but only when enough points exist (avoid degenerate models).
#
# Input (in ${THERMAL_PROJECT_ROOT}/60weeks):
#   1) 60week_NP_combined_relative_abundance_output.xlsx
#      columns: Genome, group (highT/lowT), sample columns like E01/W12/Z36...
#      NOTE: If a time point was not metagenome-sequenced, it is assumed to be 0/NA.
#   2) ENV.xlsx
#      columns: week (same naming like E01/W12/Z36...), Temperature, ...
#
# Output:
#   ${THERMAL_PROJECT_ROOT}/60weeks/out_60weeks_by_station_mgonly/
#     - station_E/, station_W/, station_Z/
#     - stability_metrics_all_stations.csv
#     - COMPARE_ratio_by_station_mgonly.png
#     - merged table: 60weeks_station_week_table_mgonly.csv
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(mgcv)
  library(segmented)
  library(stats)
  library(tibble)
})

# -----------------------------
# 0) Paths
# -----------------------------
workdir <- file.path(.release_root, "60weeks")
setwd(workdir)

abund_xlsx <- file.path(workdir, "60week_NP_combined_relative_abundance_output.xlsx")
env_xlsx   <- file.path(workdir, "ENV.xlsx")
stopifnot(file.exists(abund_xlsx), file.exists(env_xlsx))

outdir <- file.path(workdir, "out_60weeks_by_station_mgonly")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Read files
# -----------------------------
ab0 <- readxl::read_excel(abund_xlsx, sheet = 1) %>% as.data.frame()
ev0 <- readxl::read_excel(env_xlsx,   sheet = 1) %>% as.data.frame()

if (!all(c("Genome", "group") %in% names(ab0))) {
  stop("Abundance file must contain columns: Genome, group. Found: ",
       paste(names(ab0), collapse = ", "))
}
if (!all(c("week", "Temperature") %in% names(ev0))) {
  stop("ENV file must contain columns: week, Temperature. Found: ",
       paste(names(ev0), collapse = ", "))
}

# -----------------------------
# 2) Helper: parse station/week from sample name
#    Supports: E01, E1, W12, Z36, E_01, E-01 (flexible)
# -----------------------------
parse_station_week <- function(x) {
  x <- as.character(x)
  x2 <- stringr::str_replace_all(x, "[\\s_\\-]", "")
  m <- stringr::str_match(x2, "^([EWZ])0*([0-9]+)$")
  tibble::tibble(
    sample = x,
    station = m[, 2],
    week_num = suppressWarnings(as.integer(m[, 3]))
  )
}

# -----------------------------
# 3) Abundance matrix -> long -> sums by station-week-group
# -----------------------------
ab_long <- ab0 %>%
  tidyr::pivot_longer(
    cols = -c(Genome, group),
    names_to = "sample",
    values_to = "abund"
  ) %>%
  dplyr::mutate(
    group = tolower(as.character(group)),
    abund = suppressWarnings(as.numeric(abund))
  )

# Parse station/week
sample_map <- parse_station_week(unique(ab_long$sample))

ab_long <- ab_long %>%
  dplyr::left_join(sample_map, by = "sample") %>%
  dplyr::filter(!is.na(station), !is.na(week_num)) %>%
  dplyr::filter(is.finite(abund))  # drop NA

# Map group label -> HighT/LowT (adjust if needed)
ab_long <- ab_long %>%
  dplyr::mutate(group2 = dplyr::case_when(
    group %in% c("hight", "high", "warm", "warmt") ~ "HighT",
    group %in% c("lowt", "low", "cold", "coldt")  ~ "LowT",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(!is.na(group2))

# Sum across genomes
ab_sum <- ab_long %>%
  dplyr::group_by(station, week_num, group2) %>%
  dplyr::summarize(abund = sum(abund, na.rm = TRUE), .groups = "drop")

ab_week <- ab_sum %>%
  tidyr::pivot_wider(names_from = group2, values_from = abund, values_fill = 0) %>%
  dplyr::mutate(
    HighT = ifelse(is.na(HighT), 0, HighT),
    LowT  = ifelse(is.na(LowT),  0, LowT),
    AOA_total = HighT + LowT
  ) %>%
  dplyr::filter(AOA_total > 0)  # ===== KEEP ONLY metagenome points =====

# -----------------------------
# 4) ENV parse station/week the same way
# -----------------------------
ev <- ev0 %>%
  dplyr::transmute(
    sample = as.character(week),
    Temperature = suppressWarnings(as.numeric(Temperature)),
    Salinity = if ("Salinity" %in% names(ev0)) suppressWarnings(as.numeric(ev0$Salinity)) else NA_real_,
    Chla     = if ("Chla"     %in% names(ev0)) suppressWarnings(as.numeric(ev0$Chla))     else NA_real_,
    NO3      = if ("NO3"      %in% names(ev0)) suppressWarnings(as.numeric(ev0$NO3))      else NA_real_,
    NH4      = if ("NH4"      %in% names(ev0)) suppressWarnings(as.numeric(ev0$NH4))      else NA_real_
  ) %>%
  dplyr::left_join(parse_station_week(.$sample), by = "sample") %>%
  dplyr::filter(!is.na(station), !is.na(week_num)) %>%
  dplyr::arrange(station, week_num)

# -----------------------------
# 5) Merge (ab_week is backbone now, because we only analyze MG points)
# -----------------------------
dat_all <- ab_week %>%
  dplyr::left_join(ev, by = c("station", "week_num")) %>%
  dplyr::filter(is.finite(Temperature)) %>%          # must have Temperature
  dplyr::group_by(station) %>%
  dplyr::arrange(week_num, .by_group = TRUE) %>%
  dplyr::mutate(week_index = dplyr::row_number()) %>%
  dplyr::ungroup()

# epsilon per station (avoid log(0))  -- uses dplyr::select to avoid MASS conflict
eps_by_station <- dat_all %>%
  dplyr::group_by(station) %>%
  dplyr::summarize(
    min_pos = suppressWarnings(min(c(HighT[HighT > 0], LowT[LowT > 0]), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(eps = ifelse(is.finite(min_pos), min_pos/2, 1e-6)) %>%
  dplyr::select(station, eps)

dat_all <- dat_all %>%
  dplyr::left_join(eps_by_station, by = "station") %>%
  dplyr::mutate(
    ratio = log10((HighT + eps) / (LowT + eps)),
    DominantHigh = dplyr::if_else(HighT > LowT, 1L, 0L)
  )

write.csv(dat_all, file.path(outdir, "60weeks_station_week_table_mgonly.csv"), row.names = FALSE)

stations <- sort(unique(dat_all$station))

# -----------------------------
# 6) Per-station analysis (with safety guards for small n)
# -----------------------------
analyze_one_station <- function(df, station, outbase) {
  sdir <- file.path(outbase, paste0("station_", station))
  dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(df, file.path(sdir, paste0("weekly_table_", station, "_mgonly.csv")), row.names = FALSE)
  
  # Basic counts
  n_all <- nrow(df)
  n_cc  <- sum(is.finite(df$Temperature) & is.finite(df$ratio))
  n_flip_states <- length(unique(df$DominantHigh))
  
  # --- Time series plots (week_num is real time axis) ---
  p_temp <- ggplot2::ggplot(df, ggplot2::aes(x = week_num, y = Temperature)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Week number", y = "Temperature (°C)",
                  title = paste0("Station ", station, " (MG points): Temperature"))
  
  p_ratio <- ggplot2::ggplot(df, ggplot2::aes(x = week_num, y = ratio)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Week number", y = "log10((HighT+eps)/(LowT+eps))",
                  title = paste0("Station ", station, " (MG points): ratio"))
  
  ggplot2::ggsave(file.path(sdir, "TS_temperature_mgonly.png"), p_temp, width = 8, height = 4, dpi = 300)
  ggplot2::ggsave(file.path(sdir, "TS_ratio_mgonly.png"),       p_ratio, width = 8, height = 4, dpi = 300)
  
  # --- Stability metrics ---
  st <- df %>%
    dplyr::filter(is.finite(ratio)) %>%
    dplyr::summarize(
      n = dplyr::n(),
      p_flip = mean(ratio > 0, na.rm = TRUE),
      mean_ratio = mean(ratio, na.rm = TRUE),
      sd_ratio   = sd(ratio, na.rm = TRUE),
      iqr_ratio  = IQR(ratio, na.rm = TRUE),
      range_ratio = max(ratio, na.rm = TRUE) - min(ratio, na.rm = TRUE),
      cv_ratio = sd_ratio / (abs(mean_ratio) + 1e-12),
      temp_min = min(Temperature, na.rm = TRUE),
      temp_max = max(Temperature, na.rm = TRUE)
    )
  write.csv(st, file.path(sdir, "stability_metrics.csv"), row.names = FALSE)
  
  # --- GAM (only if enough points) ---
  if (n_cc >= 8) {
    # k should be small when n is small
    k_use <- ifelse(n_cc < 15, 4, 6)
    m_ratio <- mgcv::gam(ratio ~ s(Temperature, k = k_use), data = df, method = "REML")
    
    sink(file.path(sdir, "GAM_summary.txt"))
    cat("=== GAM ratio ~ s(Temperature) ===\n")
    cat("n_complete =", n_cc, "  k =", k_use, "\n\n")
    print(summary(m_ratio))
    sink()
    
    p_gam <- ggplot2::ggplot(df, ggplot2::aes(x = Temperature, y = ratio)) +
      ggplot2::geom_point() +
      ggplot2::stat_smooth(method = "gam", formula = y ~ s(x, k = k_use), se = TRUE) +
      ggplot2::labs(x = "Temperature (°C)", y = "ratio",
                    title = paste0("Station ", station, ": GAM (k=", k_use, "), MG points"))
    ggplot2::ggsave(file.path(sdir, "GAM_ratio.png"), p_gam, width = 6, height = 4.5, dpi = 300)
  } else {
    writeLines(paste0("Too few complete cases for GAM: n_complete=", n_cc, "\n"),
               con = file.path(sdir, "GAM_summary.txt"))
  }
  
  # --- CCF (only if enough sequential points; MG points may be sparse) ---
  if (n_cc >= 12) {
    cc <- stats::ccf(df$Temperature, df$ratio, lag.max = 6, plot = FALSE)
    cc_tbl <- data.frame(lag = as.integer(cc$lag), acf = as.numeric(cc$acf))
    write.csv(cc_tbl, file.path(sdir, "CCF_temp_vs_ratio.csv"), row.names = FALSE)
    
    p_ccf <- ggplot2::ggplot(cc_tbl, ggplot2::aes(x = lag, y = acf)) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_col() +
      ggplot2::labs(x = "Lag (weeks)", y = "CCF",
                    title = paste0("Station ", station, ": CCF (MG points; interpret cautiously)"))
    ggplot2::ggsave(file.path(sdir, "CCF_temp_vs_ratio.png"), p_ccf, width = 7, height = 4, dpi = 300)
  } else {
    writeLines(paste0("Too few points for CCF: n_complete=", n_cc, "\n"),
               con = file.path(sdir, "CCF_temp_vs_ratio.txt"))
  }
  
  # --- segmented breakpoint (needs enough points & variation) ---
  if (n_cc >= 15) {
    lm0 <- stats::lm(ratio ~ Temperature, data = df)
    seg_fit <- try(segmented::segmented(lm0, seg.Z = ~Temperature, psi = 18), silent = TRUE)
    
    if (!inherits(seg_fit, "try-error")) {
      psi_mat <- seg_fit$psi
      bp <- if ("Est." %in% colnames(psi_mat)) as.numeric(psi_mat[, "Est."])[1] else NA_real_
      
      sink(file.path(sdir, "SEGMENTED_summary.txt"))
      cat("=== segmented(ratio ~ Temperature) ===\n\n")
      cat("n_complete =", n_cc, "\n\n")
      cat("seg_fit$psi:\n"); print(psi_mat); cat("\n")
      cat("summary(seg_fit):\n"); print(summary(seg_fit)); cat("\n")
      cat("Estimated breakpoint (°C): ", bp, "\n")
      sink()
      
      predT <- seq(min(df$Temperature), max(df$Temperature), length.out = 200)
      predY <- predict(seg_fit, newdata = data.frame(Temperature = predT))
      
      p_seg <- ggplot2::ggplot(df, ggplot2::aes(x = Temperature, y = ratio)) +
        ggplot2::geom_point() +
        ggplot2::geom_line(data = data.frame(Temperature = predT, ratio = predY),
                           ggplot2::aes(x = Temperature, y = ratio)) +
        { if (is.finite(bp)) ggplot2::geom_vline(xintercept = bp, linetype = 2) } +
        ggplot2::labs(x = "Temperature (°C)", y = "ratio",
                      title = if (is.finite(bp)) {
                        paste0("Station ", station, ": segmented bp ≈ ", round(bp, 2), "°C (MG points)")
                      } else {
                        paste0("Station ", station, ": segmented fit (bp NA)")
                      })
      ggplot2::ggsave(file.path(sdir, "SEGMENTED_ratio_breakpoint.png"),
                      p_seg, width = 6, height = 4.5, dpi = 300)
    } else {
      writeLines("segmented failed (convergence/contrast).", con = file.path(sdir, "SEGMENTED_summary.txt"))
    }
  } else {
    writeLines(paste0("Too few points for segmented: n_complete=", n_cc, "\n"),
               con = file.path(sdir, "SEGMENTED_summary.txt"))
  }
  
  # --- logistic flip (only if both 0 and 1 exist) ---
  if (n_flip_states >= 2 && n_cc >= 10) {
    glm1 <- stats::glm(DominantHigh ~ Temperature, data = df, family = stats::binomial())
    b0 <- coef(glm1)[1]; b1 <- coef(glm1)[2]
    flipT <- as.numeric(-b0 / b1)
    
    sink(file.path(sdir, "LOGISTIC_summary.txt"))
    cat("=== Logistic: DominantHigh ~ Temperature ===\n")
    cat("n =", n_all, "  n_complete =", n_cc, "\n\n")
    print(summary(glm1))
    cat("\nFlip temperature (P=0.5): ", flipT, " °C\n")
    sink()
    
    predT <- seq(min(df$Temperature), max(df$Temperature), length.out = 200)
    predP <- predict(glm1, newdata = data.frame(Temperature = predT), type = "response")
    
    p_log <- ggplot2::ggplot(df, ggplot2::aes(x = Temperature, y = DominantHigh)) +
      ggplot2::geom_jitter(height = 0.05, width = 0, alpha = 0.6) +
      ggplot2::geom_line(data = data.frame(Temperature = predT, p = predP),
                         ggplot2::aes(x = Temperature, y = p)) +
      ggplot2::geom_vline(xintercept = flipT, linetype = 2) +
      ggplot2::labs(x = "Temperature (°C)", y = "P(HighT dominates)",
                    title = paste0("Station ", station, ": flip point ≈ ", round(flipT, 2), "°C (MG points)"))
    ggplot2::ggsave(file.path(sdir, "LOGISTIC_flip_point.png"), p_log, width = 6, height = 4.5, dpi = 300)
  } else {
    writeLines(paste0("Logistic not applicable: DominantHigh states=", n_flip_states,
                      ", n_complete=", n_cc, "\n"),
               con = file.path(sdir, "LOGISTIC_summary.txt"))
  }
  
  st %>% dplyr::mutate(station = station)
}

# -----------------------------
# 7) Run per station + summary
# -----------------------------
stability_all <- lapply(stations, function(s) {
  df_s <- dat_all %>% dplyr::filter(station == s) %>% dplyr::arrange(week_num)
  analyze_one_station(df_s, s, outdir)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::select(station, dplyr::everything())

write.csv(stability_all, file.path(outdir, "stability_metrics_all_stations.csv"), row.names = FALSE)

# Comparison plot (MG points only)
p_comp <- ggplot2::ggplot(dat_all, ggplot2::aes(x = week_num, y = ratio)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ station, ncol = 1, scales = "free_x") +
  ggplot2::labs(x = "Week number", y = "ratio",
                title = "Ratio time series by station (MG points only)")

ggplot2::ggsave(file.path(outdir, "COMPARE_ratio_by_station_mgonly.png"),
                p_comp, width = 9, height = 8, dpi = 300)

cat("\nDone!\nOutput folder: ", outdir,
    "\nStations detected (MG points only): ", paste(stations, collapse = ", "),
    "\nTotal rows analyzed: ", nrow(dat_all),
    "\n", sep = "")