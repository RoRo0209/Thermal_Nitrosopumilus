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
setwd(thermal_path("SSP585/NP_lowT"))

## =========================================================
## AOA future projection (KEEP-CHL RF ensemble) — depthsurf edition
## - Current file:  AOA_points_with_BioORACLE.tsv  (baseline predictors, surf or mean OK)
## - Future  file:  AOA_future_cov_<SCENARIO>_2080-2100_depthsurf.tsv
## - Auto-handles:
##   * baseline depthsurf vs depthmean column naming
##   * future thetao/so -> temp/sal
##   * OPTIONAL predictors missing => auto-drop (o2/no3/po4)
##   * REQUIRED predictors must exist in BOTH current & future: temp, sal, chl, si
##   * NA rows in future => predict only usable rows; keep NA for others
## Outputs:
##   out_dir/AOA_pred_future_KEEPCHL_<SCENARIO>_<PERIOD>.tsv
##   out_dir/AOA_pred_current_KEEPCHL_<SCENARIO>_<PERIOD>.tsv
## =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ranger)
  library(tibble)
  library(tidyr)
})

# -----------------------------
# Helpers
# -----------------------------
read_any_table <- function(path) {
  sz <- file.info(path)$size
  if (is.na(sz) || sz == 0) stop("EMPTY file: ", path)
  
  df <- tryCatch(read_tsv(path, show_col_types = FALSE, comment = "#", progress = FALSE),
                 error = function(e) tibble())
  if (ncol(df) == 0) {
    df <- tryCatch(read_delim(path, delim = "\t", show_col_types = FALSE, comment = "#",
                              trim_ws = TRUE, progress = FALSE),
                   error = function(e) tibble())
  }
  if (ncol(df) == 0) {
    df <- tryCatch(read_csv(path, show_col_types = FALSE, comment = "#",
                            trim_ws = TRUE, progress = FALSE),
                   error = function(e) tibble())
  }
  if (ncol(df) == 0) stop("Cannot parse (0 cols): ", path)
  df
}

detect_lonlat_cols <- function(nm) {
  lon_hits <- intersect(nm, c("Longitude","longitude","LONGITUDE","Lon","lon","LON","X","x","long","LONG"))
  lat_hits <- intersect(nm, c("Latitude","latitude","LATITUDE","Lat","lat","LAT","Y","y","lati","LATi"))
  
  lon_hits2 <- intersect(nm, c("Longitude.x","longitude.x","lon.x","Lon.x","X.x","x.x",
                               "Longitude.y","longitude.y","lon.y","Lon.y","X.y","x.y"))
  lat_hits2 <- intersect(nm, c("Latitude.x","latitude.x","lat.x","Lat.x","Y.x","y.x",
                               "Latitude.y","latitude.y","lat.y","Lat.y","Y.y","y.y"))
  
  if (length(lon_hits) > 0 && length(lat_hits) > 0) return(list(lon = lon_hits[1], lat = lat_hits[1]))
  
  lon_x <- lon_hits2[str_detect(lon_hits2, "\\.x$")]
  lat_x <- lat_hits2[str_detect(lat_hits2, "\\.x$")]
  if (length(lon_x) > 0 && length(lat_x) > 0) return(list(lon = lon_x[1], lat = lat_x[1]))
  
  if (length(lon_hits2) > 0 && length(lat_hits2) > 0) return(list(lon = lon_hits2[1], lat = lat_hits2[1]))
  
  stop("Cannot detect lon/lat cols. Available:\n", paste(nm, collapse = ", "))
}

unify_future_cov_names <- function(df) {
  nm <- names(df)
  if ("thetao" %in% nm && !"temp" %in% nm) df <- df %>% rename(temp = thetao)
  if ("so"    %in% nm && !"sal"  %in% nm) df <- df %>% rename(sal  = so)
  
  nm2 <- names(df)
  if (!"temp" %in% nm2 && any(str_detect(nm2, "^thetao_"))) {
    hit <- nm2[str_detect(nm2, "^thetao_")][1]
    df <- df %>% rename(temp = all_of(hit))
  }
  nm2 <- names(df)
  if (!"sal" %in% nm2 && any(str_detect(nm2, "^so_"))) {
    hit <- nm2[str_detect(nm2, "^so_")][1]
    df <- df %>% rename(sal = all_of(hit))
  }
  df
}

# -----------------------------
# Main function
# -----------------------------
run_projection <- function(
    base_dir,
    scenario = "SSP585",
    period = "2080-2100",
    depth_tag = "depthsurf",
    cur_file = "AOA_points_with_BioORACLE.tsv",
    future_cov_dir = "AOA_future_cov_outputs",
    out_dir = "AOA_future_projection_outputs",
    ntree = 1000,
    n_seeds = 10,
    min_node = 3
) {
  stopifnot(dir.exists(base_dir))
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(base_dir)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # -----------------------------
  # 1) Load current baseline points
  # -----------------------------
  if (!file.exists(cur_file)) stop("cur_file not found: ", file.path(base_dir, cur_file))
  cur0 <- read_tsv(cur_file, show_col_types = FALSE)
  
  y_col <- if ("AOA_log10" %in% names(cur0)) "AOA_log10" else if ("AOA_ab" %in% names(cur0)) "AOA_ab" else NA
  if (is.na(y_col)) stop("Cannot find response column AOA_log10 or AOA_ab in: ", cur_file)
  
  # baseline predictor name map: support depthsurf + depthmean + short names
  map_candidates <- list(
    temp = c("temp",
             "thetao_baseline_2000_2019_depthsurf","thetao_baseline_2000_2018_depthsurf",
             "thetao_baseline_2000_2019_depthmean","thetao_baseline_2000_2018_depthmean"),
    sal  = c("sal",
             "so_baseline_2000_2019_depthsurf","so_baseline_2000_2018_depthsurf",
             "so_baseline_2000_2019_depthmean","so_baseline_2000_2018_depthmean"),
    o2   = c("o2",
             "o2_baseline_2000_2018_depthsurf","o2_baseline_2000_2018_depthmean"),
    no3  = c("no3",
             "no3_baseline_2000_2018_depthsurf","no3_baseline_2000_2018_depthmean"),
    po4  = c("po4",
             "po4_baseline_2000_2018_depthsurf","po4_baseline_2000_2018_depthmean"),
    si   = c("si",
             "si_baseline_2000_2018_depthsurf","si_baseline_2000_2018_depthmean"),
    chl  = c("chl",
             "chl_baseline_2000_2018_depthsurf","chl_baseline_2000_2018_depthmean")
  )
  
  col_map <- sapply(names(map_candidates), function(v) {
    hits <- intersect(map_candidates[[v]], names(cur0))
    if (length(hits) == 0) NA_character_ else hits[1]
  })
  
  message("\n[CURRENT] predictor columns detected:")
  print(col_map)
  
  cur <- cur0 %>%
    transmute(
      Longitude = suppressWarnings(as.numeric(Longitude)),
      Latitude  = suppressWarnings(as.numeric(Latitude)),
      y    = suppressWarnings(as.numeric(.data[[y_col]])),
      temp = if (!is.na(col_map["temp"])) suppressWarnings(as.numeric(.data[[col_map["temp"]]])) else NA_real_,
      sal  = if (!is.na(col_map["sal"]))  suppressWarnings(as.numeric(.data[[col_map["sal"]]]))  else NA_real_,
      o2   = if (!is.na(col_map["o2"]))   suppressWarnings(as.numeric(.data[[col_map["o2"]]]))   else NA_real_,
      no3  = if (!is.na(col_map["no3"]))  suppressWarnings(as.numeric(.data[[col_map["no3"]]]))  else NA_real_,
      po4  = if (!is.na(col_map["po4"]))  suppressWarnings(as.numeric(.data[[col_map["po4"]]]))  else NA_real_,
      si   = if (!is.na(col_map["si"]))   suppressWarnings(as.numeric(.data[[col_map["si"]]]))   else NA_real_,
      chl  = if (!is.na(col_map["chl"]))  suppressWarnings(as.numeric(.data[[col_map["chl"]]]))  else NA_real_
    ) %>%
    filter(!is.na(y), !is.na(temp), !is.na(sal), !is.na(Longitude), !is.na(Latitude))
  
  message("[CURRENT] rows after filter: ", nrow(cur))
  if (nrow(cur) < 30) stop("Too few current rows after filter: ", nrow(cur))
  
  # -----------------------------
  # 2) Determine future cov file (explicit)
  # -----------------------------
  fut_file <- file.path(future_cov_dir,
                        paste0("AOA_future_cov_", scenario, "_", period, "_", depth_tag, ".tsv"))
  if (!file.exists(fut_file)) {
    stop("Future cov file not found: ", file.path(base_dir, fut_file),
         "\nExpected name: AOA_future_cov_", scenario, "_", period, "_", depth_tag, ".tsv")
  }
  message("\n[FUTURE] Using future cov file:\n  ", normalizePath(fut_file))
  
  # -----------------------------
  # 3) Load future cov and unify names
  # -----------------------------
  fut0 <- read_any_table(fut_file)
  ll <- detect_lonlat_cols(names(fut0))
  
  fut1 <- fut0 %>%
    rename(Longitude = all_of(ll$lon),
           Latitude  = all_of(ll$lat)) %>%
    unify_future_cov_names() %>%
    mutate(
      Longitude = suppressWarnings(as.numeric(Longitude)),
      Latitude  = suppressWarnings(as.numeric(Latitude))
    )
  
  cov_cols <- setdiff(names(fut1), c("Longitude","Latitude"))
  fut1[cov_cols] <- lapply(fut1[cov_cols], function(x) suppressWarnings(as.numeric(x)))
  
  fut <- fut1 %>%
    transmute(
      Longitude, Latitude,
      temp = if ("temp" %in% names(fut1)) temp else NA_real_,
      sal  = if ("sal"  %in% names(fut1)) sal  else NA_real_,
      o2   = if ("o2"   %in% names(fut1)) o2   else NA_real_,
      no3  = if ("no3"  %in% names(fut1)) no3  else NA_real_,
      po4  = if ("po4"  %in% names(fut1)) po4  else NA_real_,
      si   = if ("si"   %in% names(fut1)) si   else NA_real_,
      chl  = if ("chl"  %in% names(fut1)) chl  else NA_real_
    ) %>%
    filter(!is.na(Longitude), !is.na(Latitude))
  
  message("[FUTURE] total rows: ", nrow(fut))
  if (nrow(fut) == 0) stop("Future table has 0 rows after lon/lat parsing.")
  
  # -----------------------------
  # 4) Train + predict (KEEP CHL + SI; optional drop for others)
  # -----------------------------
  required_vars  <- c("temp","sal","chl","si")
  optional_vars  <- c("o2","no3","po4")
  
  # REQUIRED must exist in both tables (as columns)
  miss_cur_req <- setdiff(required_vars, names(cur))
  miss_fut_req <- setdiff(required_vars, names(fut))
  if (length(miss_cur_req) > 0) stop("[CURRENT] missing REQUIRED columns: ", paste(miss_cur_req, collapse=", "))
  if (length(miss_fut_req) > 0) stop("[FUTURE] missing REQUIRED columns: ", paste(miss_fut_req, collapse=", "))
  
  # For optional vars: only keep those that have at least some non-NA in BOTH
  opt_cur_ok <- optional_vars[sapply(cur[optional_vars], function(x) sum(!is.na(x)) > 0)]
  opt_fut_ok <- optional_vars[sapply(fut[optional_vars], function(x) sum(!is.na(x)) > 0)]
  opt_keep   <- intersect(opt_cur_ok, opt_fut_ok)
  
  pred_vars <- c(required_vars, opt_keep)
  
  message("\n[VARS] REQUIRED: ", paste(required_vars, collapse=", "))
  message("[VARS] OPTIONAL available (current): ", paste(opt_cur_ok, collapse=", "))
  message("[VARS] OPTIONAL available (future) : ", paste(opt_fut_ok, collapse=", "))
  message("[VARS] using predictors: ", paste(pred_vars, collapse=", "))
  
  # training complete cases
  cur_train <- cur %>% select(y, all_of(pred_vars)) %>% drop_na()
  message("[TRAIN] complete-case rows: ", nrow(cur_train))
  if (nrow(cur_train) < 30) stop("Too few complete cases for training: ", nrow(cur_train))
  
  # future usable rows (complete for pred_vars)
  ok_fut <- fut %>%
    transmute(ok = if_all(all_of(pred_vars), ~ !is.na(.))) %>%
    pull(ok)
  
  message("[FUTURE] rows usable for prediction: ", sum(ok_fut), " / ", nrow(fut))
  if (sum(ok_fut) == 0) stop("No usable future rows (all missing predictors).")
  
  mtry0 <- max(1, floor(sqrt(length(pred_vars))))
  pred_mat <- matrix(NA_real_, nrow = nrow(fut), ncol = n_seeds)
  
  for (i in seq_len(n_seeds)) {
    set.seed(i)
    fit <- ranger(
      y ~ .,
      data = cur_train,
      num.trees = ntree,
      mtry = mtry0,
      min.node.size = min_node,
      importance = "permutation"
    )
    tmp <- rep(NA_real_, nrow(fut))
    tmp[ok_fut] <- predict(fit, data = fut[ok_fut, pred_vars, drop = FALSE])$predictions
    pred_mat[, i] <- tmp
  }
  
  pred_mean <- rowMeans(pred_mat, na.rm = TRUE)
  pred_sd   <- apply(pred_mat, 1, sd, na.rm = TRUE)
  pred_cv   <- pred_sd / pmax(1e-12, abs(pred_mean))
  
  out_future <- fut %>%
    mutate(
      pred_future_mean = pred_mean,
      pred_future_sd   = pred_sd,
      pred_future_cv   = pred_cv
    )
  
  out_future_path <- file.path(out_dir, paste0("AOA_pred_future_KEEPCHL_", scenario, "_", period, ".tsv"))
  write_tsv(out_future, out_future_path)
  message("[OUT] wrote: ", normalizePath(out_future_path))
  
  # -----------------------------
  # 5) Predict CURRENT (same ensemble) for delta later
  # -----------------------------
  pred_mat_cur <- matrix(NA_real_, nrow = nrow(cur_train), ncol = n_seeds)
  
  for (i in seq_len(n_seeds)) {
    set.seed(i)
    fit <- ranger(
      y ~ .,
      data = cur_train,
      num.trees = ntree,
      mtry = mtry0,
      min.node.size = min_node
    )
    pred_mat_cur[, i] <- predict(fit, data = cur_train)$predictions
  }
  
  cur_pred_mean <- rowMeans(pred_mat_cur)
  
  cur_points <- cur %>%
    select(Longitude, Latitude, y, all_of(pred_vars)) %>%
    drop_na() %>%
    mutate(pred_current_mean = cur_pred_mean)
  
  out_cur_path <- file.path(out_dir, paste0("AOA_pred_current_KEEPCHL_", scenario, "_", period, ".tsv"))
  write_tsv(cur_points, out_cur_path)
  message("[OUT] wrote: ", normalizePath(out_cur_path))
  
  invisible(list(
    scenario = scenario,
    period = period,
    predictors_used = pred_vars,
    fut_file_used = normalizePath(fut_file),
    out_future = normalizePath(out_future_path),
    out_current = normalizePath(out_cur_path)
  ))
}

## -----------------------------
## Example usage
## -----------------------------
res_585 <- run_projection(
  base_dir  = thermal_path("SSP585/NP_lowT"),
  scenario  = "SSP585",
  period    = "2080-2100",
  depth_tag = "depthsurf"
)

print(res_585)

# res_245 <- run_projection(
#   base_dir  = thermal_path("SSP585/NP_lowT"),
#   scenario  = "SSP245",
#   period    = "2080-2100",
#   depth_tag = "depthsurf"
# )
#
# NP_highT 同理改 base_dir
# res_high_585 <- run_projection(
#   base_dir  = thermal_path("SSP585/NP_highT"),
#   scenario  = "SSP585",
#   period    = "2080-2100",
#   depth_tag = "depthsurf"
# )

