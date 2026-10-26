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
THERMAL_BIOORACLE_ROOT <- thermal_required_env("THERMAL_BIOORACLE_ROOT", directory = TRUE)
setwd(thermal_path("SSP585/NP_highT"))

## =========================================================
## 06_GLOBAL_nearshore_projection_and_maps.R  (FINAL - FIX lon/lat mutate)
## - FIX: std_lonlat_df uses .data[[]] and avoids rename(Longitude=Longitude)
## - Chunk-read TSV via base readLines + readr::read_tsv (no col_select)
## - Nearshore by dist_to_land_km <= nearshore_km AND depth >= -shelf_m
## - Plot with geom_tile (stable)
## =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ranger)
  library(ggplot2)
  library(scales)
  library(terra)
})

## -----------------------------
## 0) USER SETTINGS  (EDIT HERE)
## -----------------------------
base_dir <- thermal_path("SSP585/NP_highT")
setwd(base_dir)

scenario   <- "SSP585"
period     <- "2080-2100"
future_dir <- file.path(THERMAL_BIOORACLE_ROOT, scenario, "2080_2100_timmean")

cur_points_file <- file.path(base_dir, "AOA_points_with_BioORACLE_nearshore_filtered.tsv")  # CHANGED: use nearshore-filtered training data

global_baseline_grid_with_dist <- file.path(THERMAL_BIOORACLE_ROOT, "biooracle_nc/GLOBAL_grid_baseline_depthsurf_with_dist_km.tsv")
global_future_grid <- file.path(future_dir, paste0("GLOBAL_grid_future_", scenario, "_", period, "_depthsurf.tsv"))

out_dir <- file.path(base_dir, "AOA_future_projection_outputs", "GLOBAL_NEARSHORE_MAPS")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

nearshore_km <- 200
shelf_m      <- 200   # depth filter: keep only depth >= -shelf_m

## Bathymetry nc (auto-detect from biooracle_nc dir)
bathy_nc_dir <- file.path(THERMAL_BIOORACLE_ROOT, "biooracle_nc")
bathy_candidates <- list.files(bathy_nc_dir, pattern = "bathym|bathy|terrain", full.names = TRUE, ignore.case = TRUE)
bathy_candidates <- bathy_candidates[grepl("\\.nc$", bathy_candidates, ignore.case = TRUE)]
if (length(bathy_candidates) == 0) stop("No bathymetry .nc found in: ", bathy_nc_dir)
bathy_nc <- bathy_candidates[1]
message("[INFO] Bathymetry nc: ", basename(bathy_nc))

mode     <- "KEEPCHL"   # "KEEPCHL" or "noCHL"
nseed    <- 10
ntree    <- 1000
min_node <- 3

grid_res_deg <- 0.25
plot_res_deg <- grid_res_deg

chunk_n_lines <- 200000

## -----------------------------
## 1) Generic helpers
## -----------------------------
stopif_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  sz <- file.info(path)$size
  if (is.na(sz) || sz == 0) stop("Empty file: ", path)
  invisible(TRUE)
}

clean_names_simple <- function(nm) {
  nm2 <- trimws(nm)
  nm2 <- gsub("\u00b0", "deg", nm2)
  nm2 <- gsub("[()]", "", nm2)
  nm2 <- gsub("[^A-Za-z0-9_]+", "_", nm2)
  nm2 <- gsub("_+", "_", nm2)
  nm2 <- gsub("^_|_$", "", nm2)
  nm2
}

detect_lonlat_cols <- function(nm_orig) {
  nm_clean <- clean_names_simple(nm_orig)
  nm_low   <- tolower(nm_clean)
  
  lon_i <- which(nm_low %in% c("longitude","lon","long","x","decimal_longitude","decimallongitude","longitude_deg","lon_deg"))
  lat_i <- which(nm_low %in% c("latitude","lat","y","decimal_latitude","decimallatitude","latitude_deg","lat_deg"))
  if (length(lon_i) == 0) lon_i <- which(grepl("lon|longitud", nm_low))
  if (length(lat_i) == 0) lat_i <- which(grepl("lat|latitud", nm_low))
  if (length(lon_i) == 0 || length(lat_i) == 0) return(c(lon=NA_character_, lat=NA_character_))
  
  c(lon = nm_orig[lon_i[1]], lat = nm_orig[lat_i[1]])
}

## ✅ FIXED: robust lon/lat standardizer
std_lonlat_df <- function(df) {
  nm_orig <- names(df)
  ll <- detect_lonlat_cols(nm_orig)
  if (is.na(ll["lon"]) || is.na(ll["lat"])) {
    stop("Cannot find lon/lat columns.\nColumns:\n  ", paste(nm_orig, collapse=", "))
  }
  message("[LON/LAT] detected: lon=", ll["lon"], " | lat=", ll["lat"])
  
  # rename only if needed (avoid rename(Longitude = Longitude))
  if (ll["lon"] != "Longitude") df <- df %>% rename(Longitude = all_of(ll["lon"]))
  if (ll["lat"] != "Latitude")  df <- df %>% rename(Latitude  = all_of(ll["lat"]))
  
  # use .data[[ ]] to avoid NSE weirdness
  df %>%
    mutate(
      Longitude = suppressWarnings(as.numeric(.data[["Longitude"]])),
      Latitude  = suppressWarnings(as.numeric(.data[["Latitude"]]))
    )
}

unify_future_cov_names <- function(df) {
  nm <- names(df)
  if ("thetao" %in% nm && !"temp" %in% nm) df <- df %>% rename(temp = thetao)
  if ("so"    %in% nm && !"sal"  %in% nm) df <- df %>% rename(sal  = so)
  
  nm2 <- names(df)
  if (!"temp" %in% nm2 && any(str_detect(nm2, "^thetao_"))) {
    df <- df %>% rename(temp = all_of(nm2[str_detect(nm2, "^thetao_")][1]))
  }
  nm2 <- names(df)
  if (!"sal" %in% nm2 && any(str_detect(nm2, "^so_"))) {
    df <- df %>% rename(sal = all_of(nm2[str_detect(nm2, "^so_")][1]))
  }
  df
}

detect_dist_col <- function(nm) {
  cands <- c("dist_to_land_km","dist_km","distance_km","dist_to_coast_km",
             "dist_to_land","dist_to_coast","dist2coast_km","coast_dist_km")
  hit <- intersect(nm, cands)
  if (length(hit) > 0) return(hit[1])
  hit2 <- nm[grepl("dist", nm, ignore.case = TRUE) & grepl("km", nm, ignore.case = TRUE)]
  if (length(hit2) > 0) return(hit2[1])
  NA_character_
}

align_lonlat_to_grid <- function(df, res_deg,
                                 xlim = c(-180, 180),
                                 ylim = c(-90,  90)) {
  df %>%
    mutate(
      Longitude = floor((.data[["Longitude"]] - xlim[1]) / res_deg) * res_deg + xlim[1] + res_deg/2,
      Latitude  = floor((.data[["Latitude"]]  - ylim[1]) / res_deg) * res_deg + ylim[1] + res_deg/2
    )
}

dedup_lonlat_mean <- function(df) {
  df %>%
    group_by(Longitude, Latitude) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

tagify <- function(mode, scenario, period, nearshore_km) {
  paste0("nearshore", nearshore_km, "km_", mode, "_", scenario, "_", period)
}

## -----------------------------
## 2) RF helpers
## -----------------------------
fit_one_rf <- function(train_df, mtry, ntree, min_node, seed) {
  set.seed(seed)
  ranger(
    y ~ .,
    data = train_df,
    num.trees = ntree,
    mtry = mtry,
    min.node.size = min_node,
    splitrule = "variance"
  )
}

predict_ensemble <- function(train_df, newX, nseed = 10, ntree = 1000, min_node = 3) {
  p <- ncol(newX)
  mtry0 <- max(1, floor(sqrt(p)))
  mtry0 <- min(mtry0, p)
  pred_mat <- matrix(NA_real_, nrow = nrow(newX), ncol = nseed)
  for (i in seq_len(nseed)) {
    fit <- fit_one_rf(train_df, mtry = mtry0, ntree = ntree, min_node = min_node, seed = i)
    pred_mat[, i] <- predict(fit, data = newX)$predictions
  }
  tibble(pred_mean = rowMeans(pred_mat, na.rm = TRUE))
}

## -----------------------------
## 3) Chunked TSV reader (no col_select param)
## -----------------------------
read_tsv_selected_chunked <- function(path, keep_cols, chunk_n_lines = 200000) {
  stopif_missing(path)
  
  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  
  header <- readLines(con, n = 1)
  if (length(header) == 0) stop("Empty file header: ", path)
  
  out <- vector("list", 2000)
  oi <- 0L
  
  repeat {
    lines <- readLines(con, n = chunk_n_lines)
    if (length(lines) == 0) break
    
    txt <- c(header, lines)
    chunk <- readr::read_tsv(
      I(paste(txt, collapse = "\n")),
      show_col_types = FALSE,
      progress = FALSE
    )
    
    kc <- intersect(keep_cols, names(chunk))
    chunk <- chunk %>% select(all_of(kc))
    
    oi <- oi + 1L
    if (oi > length(out)) out <- c(out, vector("list", 2000))
    out[[oi]] <- chunk
    
    if (oi %% 10 == 0) {
      message("[CHUNK] chunks=", oi, "  rows so far=", sum(vapply(out[1:oi], nrow, integer(1))))
      gc()
    }
  }
  
  bind_rows(out[1:oi])
}

## -----------------------------
## 4) Load CURRENT points and build training table
## -----------------------------
message("[CURRENT] loading training points ...")
stopif_missing(cur_points_file)
cur_raw <- readr::read_tsv(cur_points_file, show_col_types = FALSE, progress = FALSE)
cur0 <- std_lonlat_df(cur_raw) %>% unify_future_cov_names()

y_col <- if ("AOA_log10" %in% names(cur0)) "AOA_log10" else if ("AOA_ab" %in% names(cur0)) "AOA_ab" else NA_character_
if (is.na(y_col)) stop("Cannot find response col AOA_log10 or AOA_ab in current points.")

pick1 <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}

map_candidates <- list(
  temp = c("temp","thetao_baseline_2000_2019_depthsurf","thetao_baseline_2000_2018_depthsurf"),
  sal  = c("sal","so_baseline_2000_2019_depthsurf","so_baseline_2000_2018_depthsurf"),
  chl  = c("chl","chl_baseline_2000_2018_depthsurf"),
  si   = c("si","si_baseline_2000_2018_depthsurf"),
  no3  = c("no3","no3_baseline_2000_2018_depthsurf"),
  po4  = c("po4","po4_baseline_2000_2018_depthsurf"),
  o2   = c("o2","o2_baseline_2000_2018_depthsurf")
)
col_map <- sapply(names(map_candidates), function(v) pick1(cur0, map_candidates[[v]]))
message("[CURRENT] detected predictor columns:")
print(col_map)

cur <- cur0 %>%
  transmute(
    Longitude, Latitude,
    y    = suppressWarnings(as.numeric(.data[[y_col]])),
    temp = suppressWarnings(as.numeric(.data[[col_map["temp"]]])),
    sal  = suppressWarnings(as.numeric(.data[[col_map["sal"]]])),
    chl  = if (!is.na(col_map["chl"])) suppressWarnings(as.numeric(.data[[col_map["chl"]]])) else NA_real_,
    si   = if (!is.na(col_map["si"]))  suppressWarnings(as.numeric(.data[[col_map["si"]]]))  else NA_real_,
    no3  = if (!is.na(col_map["no3"])) suppressWarnings(as.numeric(.data[[col_map["no3"]]])) else NA_real_,
    po4  = if (!is.na(col_map["po4"])) suppressWarnings(as.numeric(.data[[col_map["po4"]]])) else NA_real_,
    o2   = if (!is.na(col_map["o2"]))  suppressWarnings(as.numeric(.data[[col_map["o2"]]]))  else NA_real_
  ) %>%
  filter(is.finite(Longitude), is.finite(Latitude), is.finite(y), is.finite(temp), is.finite(sal))

if (mode == "KEEPCHL") {
  pred_vars <- c("temp","sal","chl","si","no3","po4","o2")
  required  <- c("temp","sal","chl","si")
} else {
  pred_vars <- c("temp","sal","si","no3","po4","o2")
  required  <- c("temp","sal","si")
}
pred_vars <- pred_vars[pred_vars %in% names(cur)]
if (!all(required %in% pred_vars)) stop("Missing REQUIRED predictors in current points: ", paste(required, collapse=", "))

train_df <- cur %>% select(y, all_of(pred_vars)) %>% drop_na()
message("[TRAIN] rows: ", nrow(train_df), "  predictors: ", paste(pred_vars, collapse=", "))

## -----------------------------
## 5) Load BASELINE grid WITH dist (chunked)
## -----------------------------
message("[GRID] reading BASELINE grid with dist (chunked) ...")
hdr_base <- readr::read_tsv(global_baseline_grid_with_dist, n_max = 0, show_col_types = FALSE)
base_names <- names(hdr_base)

ll_base <- detect_lonlat_cols(base_names)
if (is.na(ll_base["lon"]) || is.na(ll_base["lat"])) stop("Cannot detect lon/lat in baseline grid header.")

dist_col <- detect_dist_col(base_names)
if (is.na(dist_col)) stop("Cannot find distance column in baseline grid WITH dist.\nColumns:\n", paste(base_names, collapse=", "))
message("[DIST] baseline dist col: ", dist_col)

thetao_cols <- base_names[str_detect(base_names, "^thetao_")]
so_cols     <- base_names[str_detect(base_names, "^so_")]

need_cols_base <- unique(c(ll_base["lon"], ll_base["lat"], dist_col, pred_vars, "thetao","so", thetao_cols, so_cols))
need_cols_base <- intersect(need_cols_base, base_names)

base_raw <- read_tsv_selected_chunked(global_baseline_grid_with_dist, need_cols_base, chunk_n_lines = chunk_n_lines)

base_grid <- base_raw %>%
  std_lonlat_df() %>%
  unify_future_cov_names()

if (dist_col != "dist_to_land_km" && dist_col %in% names(base_grid)) {
  base_grid <- base_grid %>% rename(dist_to_land_km = all_of(dist_col))
}
base_grid <- base_grid %>% mutate(dist_to_land_km = suppressWarnings(as.numeric(.data[["dist_to_land_km"]])))

base_grid <- base_grid %>%
  align_lonlat_to_grid(grid_res_deg) %>%
  dedup_lonlat_mean()

message("[GRID] baseline rows after dedup: ", nrow(base_grid))

## Extract bathymetry depth for baseline grid
message("[DEPTH] extracting bathymetry from: ", basename(bathy_nc))

## Use ncdf4 to read (handles HDF5-based nc that terra/GDAL may not support)
suppressPackageStartupMessages(library(ncdf4))
nc <- ncdf4::nc_open(bathy_nc)
bathy_varname <- NULL
for (vn in c("bathymetry_mean", "bathymetry", "depth", "elevation")) {
  if (vn %in% names(nc$var)) { bathy_varname <- vn; break }
}
if (is.null(bathy_varname)) {
  # fallback: pick first non-coordinate variable
  bathy_varname <- names(nc$var)[1]
}
message("[DEPTH] using variable: ", bathy_varname)

bathy_lon <- ncdf4::ncvar_get(nc, nc$dim[[1]]$name)  # lon or longitude
bathy_lat <- ncdf4::ncvar_get(nc, nc$dim[[2]]$name)  # lat or latitude
bathy_mat <- ncdf4::ncvar_get(nc, bathy_varname)
ncdf4::nc_close(nc)

# bathy_mat is [lon, lat] or [lon, lat, time]; squeeze time if present
if (length(dim(bathy_mat)) == 3) bathy_mat <- bathy_mat[,,1]

# Build terra raster from matrix
bathy_r <- terra::rast(
  t(bathy_mat[, rev(seq_along(bathy_lat))]),
  extent = terra::ext(min(bathy_lon) - 0.025, max(bathy_lon) + 0.025,
                      min(bathy_lat) - 0.025, max(bathy_lat) + 0.025),
  crs = "EPSG:4326"
)

bathy_pts <- terra::vect(
  base_grid %>% transmute(x = Longitude, y = Latitude),
  geom = c("x", "y"), crs = "EPSG:4326"
)
base_grid$depth_m <- as.numeric(terra::extract(bathy_r, bathy_pts, method = "simple")[, 2])
rm(bathy_r, bathy_pts, bathy_mat, bathy_lon, bathy_lat); gc()

neg_frac <- mean(base_grid$depth_m < 0, na.rm = TRUE)
message("[DEPTH] negative fraction = ", round(neg_frac, 3))

base_ns <- base_grid %>%
  filter(is.finite(dist_to_land_km), dist_to_land_km <= nearshore_km)

if (neg_frac > 0.5) {
  base_ns <- base_ns %>% filter(!is.na(depth_m), depth_m < 0, depth_m >= -shelf_m)
} else {
  base_ns <- base_ns %>% filter(!is.na(depth_m), depth_m > 0, depth_m <= shelf_m)
}

message("[NEARSHORE] baseline rows (dist <= ", nearshore_km, " km & depth >= -", shelf_m, " m): ", nrow(base_ns))
if (nrow(base_ns) == 0) stop("Baseline nearshore is empty. Check dist_to_land_km / nearshore_km.")

diag_lat <- base_ns %>% count(Latitude, name="ncell") %>% arrange(desc(ncell))
message("[DIAG] top baseline nearshore lat counts:")
print(head(diag_lat, 15))

## -----------------------------
## 6) Load FUTURE grid (chunked) + JOIN dist from baseline
## -----------------------------
message("[GRID] reading FUTURE grid (chunked) ...")
hdr_fut <- readr::read_tsv(global_future_grid, n_max = 0, show_col_types = FALSE)
fut_names <- names(hdr_fut)

ll_fut <- detect_lonlat_cols(fut_names)
if (is.na(ll_fut["lon"]) || is.na(ll_fut["lat"])) stop("Cannot detect lon/lat in future grid header.")

thetao_cols2 <- fut_names[str_detect(fut_names, "^thetao_")]
so_cols2     <- fut_names[str_detect(fut_names, "^so_")]

need_cols_fut <- unique(c(ll_fut["lon"], ll_fut["lat"], pred_vars, "thetao","so", thetao_cols2, so_cols2))
need_cols_fut <- intersect(need_cols_fut, fut_names)

fut_raw <- read_tsv_selected_chunked(global_future_grid, need_cols_fut, chunk_n_lines = chunk_n_lines)

fut_grid <- fut_raw %>%
  std_lonlat_df() %>%
  unify_future_cov_names() %>%
  align_lonlat_to_grid(grid_res_deg) %>%
  dedup_lonlat_mean()

message("[GRID] future rows after dedup: ", nrow(fut_grid))

base_dist_unique <- base_grid %>% select(Longitude, Latitude, dist_to_land_km, depth_m)

n_before <- nrow(fut_grid)
fut_grid2 <- fut_grid %>% left_join(base_dist_unique, by = c("Longitude","Latitude"))
if (nrow(fut_grid2) != n_before) stop("JOIN expanded rows unexpectedly.")
fut_grid <- fut_grid2
rm(fut_grid2)

na_rate <- mean(is.na(fut_grid$dist_to_land_km))
message("[DIST] future dist NA rate after join: ", signif(na_rate, 4))
if (na_rate > 0.01) {
  stop("Too many NA dist_to_land_km after join. Grid mismatch.\nTry adjusting grid_res_deg (current=", grid_res_deg, ").")
}

fut_ns <- fut_grid %>%
  filter(is.finite(dist_to_land_km), dist_to_land_km <= nearshore_km)

if (neg_frac > 0.5) {
  fut_ns <- fut_ns %>% filter(!is.na(depth_m), depth_m < 0, depth_m >= -shelf_m)
} else {
  fut_ns <- fut_ns %>% filter(!is.na(depth_m), depth_m > 0, depth_m <= shelf_m)
}

message("[NEARSHORE] future rows (dist <= ", nearshore_km, " km & depth >= -", shelf_m, " m): ", nrow(fut_ns))
if (nrow(fut_ns) == 0) stop("Future nearshore is empty.")

## -----------------------------
## 7) Predict
## -----------------------------
base_ok <- complete.cases(base_ns[, pred_vars, drop=FALSE])
fut_ok  <- complete.cases(fut_ns[,  pred_vars, drop=FALSE])

message("[PRED] baseline usable: ", sum(base_ok), " / ", nrow(base_ns))
message("[PRED] future   usable: ", sum(fut_ok),  " / ", nrow(fut_ns))

pred_base <- rep(NA_real_, nrow(base_ns))
pred_fut  <- rep(NA_real_, nrow(fut_ns))

if (sum(base_ok) > 0) {
  pred_base[base_ok] <- predict_ensemble(train_df, base_ns[base_ok, pred_vars, drop=FALSE],
                                         nseed=nseed, ntree=ntree, min_node=min_node)$pred_mean
}
if (sum(fut_ok) > 0) {
  pred_fut[fut_ok] <- predict_ensemble(train_df, fut_ns[fut_ok, pred_vars, drop=FALSE],
                                       nseed=nseed, ntree=ntree, min_node=min_node)$pred_mean
}

current_out <- base_ns %>%
  mutate(pred_current = pred_base) %>%
  select(Longitude, Latitude, dist_to_land_km, depth_m, any_of(pred_vars), pred_current)

future_out <- fut_ns %>%
  mutate(pred_future = pred_fut) %>%
  select(Longitude, Latitude, dist_to_land_km, depth_m, any_of(pred_vars), pred_future)

delta_out <- future_out %>%
  inner_join(current_out %>% select(Longitude, Latitude, pred_current),
             by = c("Longitude","Latitude")) %>%
  mutate(
    delta_future_minus_current = pred_future - pred_current,
    delta_pct = delta_future_minus_current / pmax(1e-12, abs(pred_current))
  )

## -----------------------------
## 8) Write tables
## -----------------------------
tag <- tagify(mode, scenario, period, nearshore_km)

out_current_tsv <- file.path(out_dir, paste0("GLOBAL_nearshore_pred_current_", tag, ".tsv"))
out_future_tsv  <- file.path(out_dir, paste0("GLOBAL_nearshore_pred_future_",  tag, ".tsv"))
out_delta_tsv   <- file.path(out_dir, paste0("GLOBAL_nearshore_pred_delta_",   tag, ".tsv"))

message("[WRITE] current -> ", out_current_tsv)
write_tsv(current_out, out_current_tsv)

message("[WRITE] future  -> ", out_future_tsv)
write_tsv(future_out,  out_future_tsv)

message("[WRITE] delta   -> ", out_delta_tsv)
write_tsv(delta_out,   out_delta_tsv)

message("[OUT] wrote tables ✅")

## -----------------------------
## 9) Plot maps (PDF) with geom_tile
## -----------------------------
theme_pub_map <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(size = base_size),
      legend.text  = element_text(size = base_size - 1),
      plot.title = element_text(size = base_size + 1, face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

world_df <- ggplot2::map_data("world")

bin_for_plot <- function(df, zcol, res_deg) {
  df %>%
    select(Longitude, Latitude, !!sym(zcol)) %>%
    filter(is.finite(Longitude), is.finite(Latitude), is.finite(.data[[zcol]])) %>%
    rename(z = !!sym(zcol)) %>%
    align_lonlat_to_grid(res_deg) %>%
    group_by(Longitude, Latitude) %>%
    summarise(z = mean(z, na.rm = TRUE), .groups = "drop")
}

plot_tile_world_pdf <- function(df, zcol, title, out_pdf,
                                res_deg = 0.25,
                                clip_quant = c(0.02, 0.98),
                                delta_mode = FALSE,
                                base_size = 11) {
  
  d <- bin_for_plot(df, zcol, res_deg)
  
  if (delta_mode) {
    q <- quantile(abs(d$z), 0.98, na.rm = TRUE)
    lims <- c(-q, q)
  } else {
    q <- quantile(d$z, probs = clip_quant, na.rm = TRUE)
    lims <- c(q[1], q[2])
  }
  
  p <- ggplot() +
    geom_rect(aes(xmin=-180, xmax=180, ymin=-90, ymax=90), fill="white", color=NA) +
    geom_tile(data = d, aes(x = Longitude, y = Latitude, fill = z),
              width = res_deg, height = res_deg) +
    geom_polygon(data = world_df, aes(x = long, y = lat, group = group),
                 fill = "grey95", color = "grey70", linewidth = 0.15) +
    coord_quickmap(xlim = c(-180, 180), ylim = c(-70, 80), expand = FALSE) +
    labs(title = title, fill = zcol) +
    theme_pub_map(base_size = base_size)
  
  if (delta_mode) {
    p <- p + scale_fill_gradient2(limits = lims, oob = squish, midpoint = 0)
  } else {
    p <- p + scale_fill_viridis_c(limits = lims, oob = squish)
  }
  
  if (capabilities("cairo")) {
    ggsave(out_pdf, p, width = 11, height = 5.5, device = cairo_pdf, bg = "white")
  } else {
    ggsave(out_pdf, p, width = 11, height = 5.5, bg = "white")
  }
  p
}

map_current_pdf <- file.path(out_dir, paste0("Map_nearshore_current_", tag, ".pdf"))
map_future_pdf  <- file.path(out_dir, paste0("Map_nearshore_future_",  tag, ".pdf"))
map_delta_pdf   <- file.path(out_dir, paste0("Map_nearshore_delta_",   tag, ".pdf"))

plot_tile_world_pdf(current_out, "pred_current",
                    paste0("GLOBAL Nearshore AOA (CURRENT)  ", mode, " | <= ", nearshore_km, " km"),
                    map_current_pdf, res_deg = plot_res_deg, delta_mode = FALSE)

plot_tile_world_pdf(future_out, "pred_future",
                    paste0("GLOBAL Nearshore AOA (FUTURE)  ", scenario, " ", period, "  ", mode, " | <= ", nearshore_km, " km"),
                    map_future_pdf, res_deg = plot_res_deg, delta_mode = FALSE)

plot_tile_world_pdf(delta_out, "delta_future_minus_current",
                    paste0("GLOBAL Nearshore AOA (DELTA future-current)  ", scenario, " ", period, "  ", mode, " | <= ", nearshore_km, " km"),
                    map_delta_pdf, res_deg = plot_res_deg, delta_mode = TRUE)

message("\nDONE ✅ outputs:")
message("  - ", normalizePath(out_current_tsv))
message("  - ", normalizePath(out_future_tsv))
message("  - ", normalizePath(out_delta_tsv))
message("  - ", normalizePath(map_current_pdf))
message("  - ", normalizePath(map_future_pdf))
message("  - ", normalizePath(map_delta_pdf))
