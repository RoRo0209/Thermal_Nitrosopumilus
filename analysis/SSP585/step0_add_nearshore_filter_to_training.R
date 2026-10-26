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
setwd(thermal_path("SSP585"))

## =========================================================
## step0_add_nearshore_filter_to_training.R
## Add depth_m and dist_to_land_km to training data, then filter:
##   depth >= -200 m  AND  dist_to_land_km <= 200 km
## Outputs: AOA_points_with_BioORACLE_nearshore_filtered.tsv
## =========================================================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ncdf4)
  library(terra)
  library(data.table)
})

## --- CONFIG ---
datasets <- list(
  list(
    name     = "NP_lowT",
    base_dir = thermal_path("SSP585/NP_lowT"),
    in_file  = "AOA_points_with_BioORACLE.tsv"
  ),
  list(
    name     = "NP_highT",
    base_dir = thermal_path("SSP585/NP_highT"),
    in_file  = "AOA_points_with_BioORACLE.tsv"
  )
)

bathy_nc_path <- file.path(THERMAL_BIOORACLE_ROOT, "biooracle_nc/terrain_characteristics.nc")
global_dist_file <- file.path(THERMAL_BIOORACLE_ROOT, "biooracle_nc/GLOBAL_grid_baseline_depthsurf_with_dist_km.tsv")

nearshore_km <- 200
shelf_m      <- 200

## --- 1) Load bathymetry raster from nc (same as step8) ---
message("[BATHY] Reading: ", basename(bathy_nc_path))
nc <- nc_open(bathy_nc_path)
bathy_varname <- NULL
for (vn in c("bathymetry_mean", "bathymetry", "depth", "elevation")) {
  if (vn %in% names(nc$var)) { bathy_varname <- vn; break }
}
if (is.null(bathy_varname)) bathy_varname <- names(nc$var)[1]
message("[BATHY] variable: ", bathy_varname)

bathy_lon <- ncvar_get(nc, nc$dim[[1]]$name)
bathy_lat <- ncvar_get(nc, nc$dim[[2]]$name)
bathy_mat <- ncvar_get(nc, bathy_varname)
nc_close(nc)

if (length(dim(bathy_mat)) == 3) bathy_mat <- bathy_mat[,,1]

bathy_r <- rast(
  t(bathy_mat[, rev(seq_along(bathy_lat))]),
  extent = ext(min(bathy_lon) - 0.025, max(bathy_lon) + 0.025,
               min(bathy_lat) - 0.025, max(bathy_lat) + 0.025),
  crs = "EPSG:4326"
)
rm(bathy_mat, bathy_lon, bathy_lat); gc()

## --- 2) Load dist_to_land_km grid for nearest-neighbor lookup ---
message("[DIST] Reading dist_to_land_km from global grid ...")
dist_dt <- fread(global_dist_file,
                 select = c("Longitude", "Latitude", "dist_to_land_km"),
                 showProgress = FALSE)
setkey(dist_dt, Longitude, Latitude)
message("[DIST] Global grid rows: ", nrow(dist_dt))

## Nearest-neighbor function: for each query point, find closest grid cell
nn_dist_lookup <- function(query_lon, query_lat, grid_dt, grid_res = 0.05) {
  ## For each point, search within a small neighborhood
  n <- length(query_lon)
  result <- numeric(n)
  for (i in seq_len(n)) {
    ql <- query_lon[i]
    qa <- query_lat[i]
    ## Search within +/- 2 grid cells
    tol <- grid_res * 2
    sub <- grid_dt[Longitude >= (ql - tol) & Longitude <= (ql + tol) &
                   Latitude >= (qa - tol) & Latitude <= (qa + tol)]
    if (nrow(sub) == 0) {
      ## Widen search
      tol <- grid_res * 5
      sub <- grid_dt[Longitude >= (ql - tol) & Longitude <= (ql + tol) &
                     Latitude >= (qa - tol) & Latitude <= (qa + tol)]
    }
    if (nrow(sub) == 0) {
      result[i] <- NA_real_
    } else {
      d2 <- (sub$Longitude - ql)^2 + (sub$Latitude - qa)^2
      result[i] <- sub$dist_to_land_km[which.min(d2)]
    }
  }
  result
}

## --- 3) Process each dataset ---
for (ds in datasets) {
  message("\n========== ", ds$name, " ==========")
  in_path <- file.path(ds$base_dir, ds$in_file)
  if (!file.exists(in_path)) stop("Required training input missing: ", in_path)

  dat <- read_tsv(in_path, show_col_types = FALSE)
  message("Input rows: ", nrow(dat))

  ## 3a) Extract bathymetry
  pts <- vect(dat, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
  dat$depth_m <- as.numeric(extract(bathy_r, pts, method = "simple")[, 2])

  ## 3b) Extract dist_to_land_km via nearest-neighbor
  dat$dist_to_land_km <- nn_dist_lookup(dat$Longitude, dat$Latitude, dist_dt)

  message("depth_m range: ", paste(round(range(dat$depth_m, na.rm=TRUE), 1), collapse=" ~ "))
  message("depth_m NA: ", sum(is.na(dat$depth_m)), " / ", nrow(dat))
  message("dist_to_land_km range: ", paste(round(range(dat$dist_to_land_km, na.rm=TRUE), 1), collapse=" ~ "))
  message("dist_to_land_km NA: ", sum(is.na(dat$dist_to_land_km)), " / ", nrow(dat))

  ## Print diagnostic before filtering
  message("\n--- Filter diagnostic ---")
  message("depth_m >= -200: ", sum(dat$depth_m >= -shelf_m, na.rm=TRUE))
  message("dist_to_land_km <= 200: ", sum(dat$dist_to_land_km <= nearshore_km, na.rm=TRUE))
  message("Both conditions: ", sum(dat$depth_m >= -shelf_m & dat$dist_to_land_km <= nearshore_km, na.rm=TRUE))

  ## 3c) Apply nearshore filter
  neg_frac <- mean(dat$depth_m < 0, na.rm = TRUE)
  message("Negative depth fraction: ", round(neg_frac, 3))

  if (neg_frac > 0.5) {
    dat_ns <- dat %>%
      filter(!is.na(depth_m), depth_m <= 0, depth_m >= -shelf_m,
             !is.na(dist_to_land_km), dist_to_land_km <= nearshore_km)
  } else {
    dat_ns <- dat %>%
      filter(!is.na(depth_m), depth_m >= 0, depth_m <= shelf_m,
             !is.na(dist_to_land_km), dist_to_land_km <= nearshore_km)
  }

  message("\nAfter nearshore filter (dist<=", nearshore_km, "km & depth>=-", shelf_m, "m): ",
          nrow(dat_ns), " / ", nrow(dat), " rows")

  ## 3d) Save filtered + also save full with depth/dist info for inspection
  out_path <- file.path(ds$base_dir, "AOA_points_with_BioORACLE_nearshore_filtered.tsv")
  write_tsv(dat_ns, out_path)
  message("Wrote: ", out_path)

  out_full <- file.path(ds$base_dir, "AOA_points_with_BioORACLE_with_depth_dist.tsv")
  write_tsv(dat, out_full)
  message("Wrote (full with depth/dist): ", out_full)
}

message("\nDONE")
rm(bathy_r, dist_dt); gc()
