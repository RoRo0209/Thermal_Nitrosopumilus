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
setwd(thermal_path("SSP585/NP_lowT"))

## =========================================================
## AOA future cov extraction from LOCAL Bio-ORACLE NetCDF (CDO-processed)
## - Input : NP_lowT_sample_abundance.tsv (Longitude, Latitude)
## - Read  : LOCAL .nc already processed to 2080-2100 timmean (single layer)
## - Output: AOA_future_cov_<SCENARIO>_2080-2100_depthsurf.tsv
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(terra)
})

setwd(thermal_path("SSP585/NP_lowT"))

## -----------------------------
## 0) Paths / settings (EDIT)
## -----------------------------
in_points <- "NP_lowT_sample_abundance.tsv"

SCENARIO  <- "SSP585"         # <- 改成 "SSP245" 或 "SSP585"
PERIOD    <- "2080-2100"
DEPTH_TAG <- "depthsurf"

# 关键：指向 CDO 输出的文件夹（里面是 *_2080_2100_timmean.nc）
scenario_dir <- file.path(THERMAL_BIOORACLE_ROOT, SCENARIO, "2080_2100_timmean")
stopifnot(dir.exists(scenario_dir))

out_dir  <- "AOA_future_cov_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, paste0(
  "AOA_future_cov_", SCENARIO, "_", PERIOD, "_", DEPTH_TAG, ".tsv"
))

extract_method <- "near"  # "near" or "bilinear"

## -----------------------------
## 1) Read points
## -----------------------------
pts_df <- read_tsv(in_points, show_col_types = FALSE) %>%
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude  = as.numeric(Latitude)
  )

stopifnot(all(c("Longitude","Latitude") %in% names(pts_df)))

pts <- pts_df %>%
  transmute(.id = row_number(), x = Longitude, y = Latitude) %>%
  terra::vect(geom = c("x","y"), crs = "EPSG:4326")

## -----------------------------
## 2) List nc files + match variables
## -----------------------------
nc_files <- list.files(scenario_dir, pattern = "\\.nc$", full.names = TRUE, ignore.case = TRUE)
if (length(nc_files) == 0) stop("No .nc files found in: ", scenario_dir)

# 输出列名（建议与你后续未来预测脚本一致）
want_vars <- c("thetao","so","o2","no3","po4","si","chl")

file_patterns <- list(
  thetao = "thetao|oceantemperature|temperature|temp",
  so     = "(^|[^a-z0-9])so([^a-z0-9]|$)|salinity|sal",
  o2     = "(^|[^a-z0-9])o2([^a-z0-9]|$)|oxygen",
  no3    = "(^|[^a-z0-9])no3([^a-z0-9]|$)|nitrate",
  po4    = "(^|[^a-z0-9])po4([^a-z0-9]|$)|phosphate",
  si     = "(^|[^a-z0-9])si([^a-z0-9]|$)|silicate",
  chl    = "(^|[^a-z0-9])chl([^a-z0-9]|$)|chlorophyll"
)

pick_nc_for_var <- function(var) {
  pat <- file_patterns[[var]]
  hit <- nc_files[str_detect(tolower(basename(nc_files)), pat)]
  if (length(hit) == 0) return(NA_character_)
  
  # 关键：优先选 timmean + 2080/2090/2100 + depthsurf
  hit2 <- hit[order(
    str_detect(tolower(basename(hit)), "timmean"),
    str_detect(tolower(basename(hit)), "2080|2090|2100|2080_2100"),
    str_detect(tolower(basename(hit)), tolower(DEPTH_TAG)),
    decreasing = TRUE
  )]
  
  hit2[1]
}

var2file <- setNames(vapply(want_vars, pick_nc_for_var, character(1)), want_vars)

message("[INFO] scenario_dir: ", scenario_dir)
message("[INFO] Matched NC files:")
print(basename(var2file))

if (any(is.na(var2file))) {
  missing <- names(var2file)[is.na(var2file)]
  stop(
    "Cannot find .nc for: ", paste(missing, collapse = ", "),
    "\nPlease check filenames in: ", scenario_dir,
    "\nTip: run `basename(list.files(scenario_dir, pattern='\\\\.nc$', full.names=TRUE))` and adjust file_patterns."
  )
}

## -----------------------------
## 3) Helpers (read layer + lon fix + extract)
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
  
  if (lon_min >= 0 && lon_max > 180 && any(pts_lon < 0, na.rm = TRUE)) {
    pts_xy2 <- pts_xy
    pts_xy2[,1] <- ifelse(pts_xy2[,1] < 0, pts_xy2[,1] + 360, pts_xy2[,1])
    pts2 <- terra::vect(pts_xy2, crs = terra::crs(pts_vect))
    return(list(r = r, pts = pts2, mode = "points shifted to 0-360"))
  }

  if (lon_min < 0 && lon_max <= 180 && any(pts_lon > 180, na.rm = TRUE)) {
    pts_xy2 <- pts_xy
    pts_xy2[,1] <- ifelse(pts_xy2[,1] > 180, pts_xy2[,1] - 360, pts_xy2[,1])
    pts2 <- terra::vect(pts_xy2, crs = terra::crs(pts_vect))
    return(list(r = r, pts = pts2, mode = "points shifted to -180~180"))
  }
  
  list(r = r, pts = pts_vect, mode = "no shift")
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
  
  list(values = as.numeric(v), lon_mode = fx$mode, extent = terra::ext(r2))
}

## -----------------------------
## 4) Extract
## -----------------------------
out_df <- pts_df %>% select(Longitude, Latitude)

for (vname in want_vars) {
  f <- var2file[[vname]]
  message("\n[EXTRACT] ", vname, " from: ", basename(f))
  
  res <- extract_one(f, pts, method = extract_method)
  vals <- res$values
  
  message("  lon handling: ", res$lon_mode,
          " | raster lon extent=", round(res$extent[1],2), "~", round(res$extent[2],2))
  message("  NA rate: ", round(mean(is.na(vals)), 3),
          " | range: ", signif(suppressWarnings(min(vals, na.rm = TRUE)), 4),
          " ~ ", signif(suppressWarnings(max(vals, na.rm = TRUE)), 4))
  
  out_df[[vname]] <- vals
}

## -----------------------------
## 5) Write
## -----------------------------
write_tsv(out_df, out_file)
message("\nDone ✅ Wrote:\n  ", out_file)
message("Columns: ", paste(names(out_df), collapse=", "))
