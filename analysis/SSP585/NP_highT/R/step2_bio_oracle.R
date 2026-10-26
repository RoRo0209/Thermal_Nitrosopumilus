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
## step2_local_BioORACLE_extract_from_NC_NEARSHORE.R
## - LOCAL .nc extraction
## - FIX: filter nearshore by bathymetry depth <= shelf_m (200 m) in step2
## Output: AOA_points_with_BioORACLE_nearshore.tsv
## =========================================================

## =========================================================
## step2_local_BioORACLE_extract_NEARSHORE_if_possible.R
## - Extract Bio-ORACLE predictors from LOCAL .nc files (terra)
## - If bathymetry (bathy/bathym) .nc exists: filter nearshore depth <= shelf_m
## - If NOT exists: skip nearshore filter (still extract predictors)
##
## Input : NP_lowT_sample_abundance.tsv
## Output:
##   1) AOA_points_with_BioORACLE.tsv                  (all points with predictors)
##   2) AOA_points_with_BioORACLE_nearshore.tsv        (if bathy exists; else same as all but with warning)
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(terra)
})

setwd(thermal_path("SSP585/NP_highT"))

in_points <- "NP_highT_sample_abundance.tsv"

out_all  <- "AOA_points_with_BioORACLE.tsv"
out_near <- "AOA_points_with_BioORACLE_nearshore.tsv"

## -----------------------------
## 0) SETTINGS
## -----------------------------
nc_dir   <- file.path(THERMAL_BIOORACLE_ROOT, "biooracle_nc")  # <- 改成你的真实目录（你现在这个是对的）
shelf_m  <- 200
extract_method <- "near"  # "near" or "bilinear"

if (!dir.exists(nc_dir)) stop("nc_dir does not exist: ", nc_dir)

## -----------------------------
## 1) Read points
## -----------------------------
pts_df <- read_tsv(in_points, show_col_types = FALSE) %>%
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude  = as.numeric(Latitude)
  )

stopifnot(all(c("Longitude","Latitude") %in% names(pts_df)))

make_vect <- function(df) {
  df %>%
    transmute(.id = row_number(), x = Longitude, y = Latitude) %>%
    terra::vect(geom = c("x","y"), crs = "EPSG:4326")
}

pts <- make_vect(pts_df)

## -----------------------------
## 2) List local nc files
## -----------------------------
nc_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE, ignore.case = TRUE)
if (length(nc_files) == 0) stop("No .nc files found in: ", nc_dir)

## predictors you want
want_vars <- c("chl","o2","no3","po4","si","temp","sal")

file_patterns <- list(
  chl  = "chl|chlorophyll",
  o2   = "o2|oxygen",
  no3  = "no3|nitrate",
  po4  = "po4|phosphate",
  si   = "si|silicate",
  temp = "thetao|temp|temperature|ocean_temperature",
  sal  = "so|sal|salinity"
)

pick_nc_for_var <- function(var, nc_files) {
  pat <- file_patterns[[var]]
  hit <- nc_files[str_detect(tolower(basename(nc_files)), pat)]
  if (length(hit) == 0) return(NA_character_)
  
  # prefer baseline + depthsurf/depthmean (any climatology layer)
  hit2 <- hit[order(
    str_detect(tolower(basename(hit)), "baseline"),
    str_detect(tolower(basename(hit)), "depthmean|depth_mean|depth-mean|depthsurf|depth_surface|surface|surf"),
    decreasing = TRUE
  )]
  hit2[1]
}

var2file <- setNames(vapply(want_vars, pick_nc_for_var, character(1), nc_files = nc_files), want_vars)

message("[INFO] Predictor files matched:")
print(var2file)

if (any(is.na(var2file))) {
  missing <- names(var2file)[is.na(var2file)]
  stop(
    "Cannot find local .nc for: ", paste(missing, collapse = ", "),
    "\nPlease check filenames under: ", nc_dir
  )
}

## -----------------------------
## 3) Bathymetry file (optional)
## IMPORTANT: DO NOT use 'depth' keyword because your files contain 'depthsurf'
## -----------------------------
bathy_hit <- nc_files[str_detect(tolower(basename(nc_files)), "bathym|bathy|terrain")]

bathy_file <- NA_character_
if (length(bathy_hit) > 0) {
  # if multiple, just take first (you can refine later)
  bathy_file <- bathy_hit[1]
  message("[INFO] Bathymetry file found: ", basename(bathy_file))
} else {
  message("[WARN] No bathymetry/bathy .nc found in nc_dir; nearshore filtering will be SKIPPED.")
  message("[WARN] If you want depth<=200m nearshore, download Bio-ORACLE bathymetry layer (.nc) and put it in: ", nc_dir)
}

## dist_coast file (optional but recommended for nearshore filter)
coast_hit <- nc_files[str_detect(tolower(basename(nc_files)), "dist_coast|distcoast|coast_dist")]
coast_file <- NA_character_
if (length(coast_hit) > 0) {
  coast_file <- coast_hit[1]
  message("[INFO] Distance-to-coast file found: ", basename(coast_file))
} else {
  message("[WARN] No dist_coast .nc found in nc_dir; dist_coast filter will be SKIPPED.")
}

## -----------------------------
## 4) Helpers
## -----------------------------
read_first_layer <- function(nc_path) {
  r <- terra::rast(nc_path)
  if (terra::nlyr(r) > 1) {
    nms <- names(r)
    idx <- which(str_detect(tolower(nms), "mean"))
    if (length(idx) > 0) r <- r[[idx[1]]] else r <- r[[1]]
  }
  r
}

fix_lon_if_needed <- function(r, pts_vect) {
  ex <- terra::ext(r)
  lon_min <- ex[1]; lon_max <- ex[2]
  pts_xy <- terra::crds(pts_vect)
  pts_lon <- pts_xy[,1]
  
  # raster 0-360, points -180~180
  if (lon_min >= 0 && lon_max > 180 && any(pts_lon < 0, na.rm = TRUE)) {
    pts_xy2 <- pts_xy
    pts_xy2[,1] <- ifelse(pts_xy2[,1] < 0, pts_xy2[,1] + 360, pts_xy2[,1])
    pts2 <- terra::vect(pts_xy2, crs = terra::crs(pts_vect))
    return(list(r = r, pts = pts2))
  }

  # raster -180~180, points 0-360
  if (lon_min < 0 && lon_max <= 180 && any(pts_lon > 180, na.rm = TRUE)) {
    pts_xy2 <- pts_xy
    pts_xy2[,1] <- ifelse(pts_xy2[,1] > 180, pts_xy2[,1] - 360, pts_xy2[,1])
    pts2 <- terra::vect(pts_xy2, crs = terra::crs(pts_vect))
    return(list(r = r, pts = pts2))
  }
  
  list(r = r, pts = pts_vect)
}

extract_one <- function(nc_path, pts_vect, method = c("near","bilinear")) {
  method <- match.arg(method)
  r <- read_first_layer(nc_path)
  
  if (is.na(terra::crs(r)) || terra::crs(r) == "") {
    suppressWarnings(terra::crs(r) <- "EPSG:4326")
  }
  
  fx <- fix_lon_if_needed(r, pts_vect)
  r2 <- fx$r
  pts2 <- fx$pts
  
  if (method == "near") {
    v <- terra::extract(r2, pts2, method = "simple")[,2]
  } else {
    v <- terra::extract(r2, pts2, method = "bilinear")[,2]
  }
  as.numeric(v)
}

## -----------------------------
## 5) Extract predictors for ALL points
## -----------------------------
out_df_all <- pts_df

for (vname in want_vars) {
  message("\n[EXTRACT] ", vname, " from: ", basename(var2file[[vname]]))
  vals <- extract_one(var2file[[vname]], pts, method = extract_method)
  out_df_all[[vname]] <- vals
  message("  NA rate: ", round(mean(is.na(vals)), 3),
          " | range: ", signif(suppressWarnings(min(vals, na.rm = TRUE)), 4),
          " ~ ", signif(suppressWarnings(max(vals, na.rm = TRUE)), 4))
}

write_tsv(out_df_all, out_all)
message("\n[DONE] Wrote ALL points: ", out_all, " (n=", nrow(out_df_all), ")")

## -----------------------------
## 6) Nearshore filter if bathy exists
## -----------------------------
out_df_near <- out_df_all

if (!is.na(bathy_file) || !is.na(coast_file)) {
  tmp <- out_df_all

  ## --- depth filter ---
  if (!is.na(bathy_file)) {
    message("\n[STEP] Extract bathymetry and filter depth >= -", shelf_m, " m ...")
    depth_raw <- extract_one(bathy_file, pts, method = "near")
    tmp <- tmp %>% mutate(depth_m = depth_raw)

    neg_frac <- mean(tmp$depth_m < 0, na.rm = TRUE)
    message("[INFO] depth negative fraction = ", round(neg_frac, 3))

    if (neg_frac > 0.5) {
      tmp <- tmp %>% filter(!is.na(depth_m), depth_m < 0, depth_m >= -shelf_m)
    } else {
      tmp <- tmp %>% filter(!is.na(depth_m), depth_m > 0, depth_m <= shelf_m)
    }
    message("[INFO] After depth filter: ", nrow(tmp), " points")
  }

  ## --- dist_coast filter ---
  if (!is.na(coast_file)) {
    message("[STEP] Extract dist_coast and filter <= 200 km ...")
    dist_raw <- extract_one(coast_file, make_vect(tmp), method = "near")
    tmp <- tmp %>% mutate(dist_coast_km = dist_raw)

    ## Bio-ORACLE dist_coast is in metres; convert if needed
    if (median(tmp$dist_coast_km, na.rm = TRUE) > 1000) {
      message("[INFO] dist_coast appears in metres, converting to km")
      tmp <- tmp %>% mutate(dist_coast_km = dist_coast_km / 1000)
    }

    tmp <- tmp %>% filter(!is.na(dist_coast_km), dist_coast_km <= 200)
    message("[INFO] After dist_coast <= 200 km filter: ", nrow(tmp), " points")
  }

  if (nrow(tmp) == 0) stop("No nearshore points left after filtering; check bathy/dist_coast files and thresholds.")

  out_df_near <- tmp
} else {
  message("\n[WARN] Nearshore filtering skipped (no bathymetry or dist_coast file). Writing out_near same as out_all.")
}

write_tsv(out_df_near, out_near)
message("[DONE] Wrote NEARSHORE: ", out_near, " (n=", nrow(out_df_near), ")")

## -----------------------------
## 7) QC summary
## -----------------------------
pred_cols <- want_vars
na_rate_all <- sapply(out_df_all[pred_cols], function(x) mean(is.na(x)))
message("\nNA rates (ALL points predictors):")
print(round(na_rate_all, 3))

na_rate_near <- sapply(out_df_near[pred_cols], function(x) mean(is.na(x)))
message("\nNA rates (NEARSHORE predictors):")
print(round(na_rate_near, 3))

message("\nALL DONE ✅")
