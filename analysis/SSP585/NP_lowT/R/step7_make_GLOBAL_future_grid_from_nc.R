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
## 03_make_GLOBAL_future_grid_from_nc.R  (FIXED + SAFE)
## - Flatten Bio-ORACLE future timmean nc into a global OCEAN-only grid table
## - Uses ONE template raster to define:
##   * grid geometry (ncol/nrow/res/extent)
##   * ocean cells (non-NA in template)
## - Extracts values ONLY on those ocean cells -> no length mismatch
##
## Input dir:
##
## Output:
##   GLOBAL_grid_future_SSP585_2080-2100_depthsurf.tsv
## =========================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

## -----------------------------
## 0) EDIT (already set for you)
## -----------------------------
scenario_dir <- file.path(THERMAL_BIOORACLE_ROOT, "SSP585/2080_2100_timmean")
scenario <- "SSP585"
period   <- "2080-2100"
depth_tag <- "depthsurf"

out_tsv <- file.path(scenario_dir,
                     paste0("GLOBAL_grid_future_", scenario, "_", period, "_", depth_tag, ".tsv"))

## -----------------------------
## 1) Choose a template nc (must be on the target grid)
## -----------------------------
# 优先找 thetao 的 depthsurf timmean 文件（你之前生成过类似：thetao_ssp585_2020_2100_depthsurf_2080_2100_timmean.nc）
candidates <- list.files(scenario_dir, pattern = "thetao.*depthsurf.*timmean.*\\.nc$", full.names = TRUE, ignore.case = TRUE)

if (length(candidates) == 0) {
  # 兜底：随便找一个 nc 做模板（必须是同网格）
  candidates <- list.files(scenario_dir, pattern = ".*depthsurf.*timmean.*\\.nc$", full.names = TRUE, ignore.case = TRUE)
}
if (length(candidates) == 0) {
  stop("No template nc found under: ", scenario_dir,
       "\nNeed something like: thetao_*depthsurf*_timmean.nc")
}

template_nc <- candidates[1]
message("[TEMPLATE] ", basename(template_nc))

r0 <- rast(template_nc)
if (nlyr(r0) > 1) r0 <- r0[[1]]
crs(r0) <- "EPSG:4326"

ncells <- ncell(r0)
message("[GRID] ncol=", ncol(r0), " nrow=", nrow(r0), " ncell=", ncells)

## -----------------------------
## 2) Define OCEAN cells as non-NA in template
## -----------------------------
v0 <- values(r0, mat = FALSE)
ocean_cells <- which(!is.na(v0))

message("[OCEAN] ocean cells (non-NA in template) = ", length(ocean_cells))

if (length(ocean_cells) < 1000000) {
  message("[WARN] Ocean cells seem unusually low; check if template contains too many NA.")
}

xy <- as.data.frame(xyFromCell(r0, ocean_cells))
colnames(xy) <- c("Longitude", "Latitude")

grid <- as_tibble(xy)

rm(v0, xy); gc()

## -----------------------------
## 3) Helper: read a variable nc and extract values on ocean cells
## -----------------------------
extract_on_ocean <- function(nc_file, template_r, ocean_cells, method = "bilinear") {
  r <- rast(nc_file)
  if (nlyr(r) > 1) r <- r[[1]]
  crs(r) <- "EPSG:4326"
  
  # If geometry mismatch, resample to template
  same_geom <- (ncol(r) == ncol(template_r) &&
                  nrow(r) == nrow(template_r) &&
                  all.equal(ext(r), ext(template_r)) == TRUE)
  
  if (!same_geom) {
    message("[RESAMPLE] ", basename(nc_file), " -> template grid (", method, ")")
    r <- resample(r, template_r, method = method)
  }
  
  vals <- values(r, mat = FALSE)[ocean_cells]
  vals
}

## -----------------------------
## 4) Map variables -> file patterns (auto-find)
## -----------------------------
# 你这个目录里通常是：thetao / so / chl / no3 / po4 / si / o2 的 timmean nc
# 下面用 pattern 自动匹配，不要写死文件名（更稳）
file_patterns <- list(
  temp = "thetao.*timmean.*\\.nc$",
  sal  = "(^|/)so.*timmean.*\\.nc$",
  chl  = "chl.*timmean.*\\.nc$",
  no3  = "no3.*timmean.*\\.nc$",
  po4  = "po4.*timmean.*\\.nc$",
  si   = "(^|/)si.*timmean.*\\.nc$",
  o2   = "(^|/)o2.*timmean.*\\.nc$"
)

# list all nc files
all_nc <- list.files(scenario_dir, pattern = "\\.nc$", full.names = TRUE, ignore.case = TRUE)
if (length(all_nc) == 0) stop("No .nc files under: ", scenario_dir)

pick_file <- function(pattern) {
  hit <- all_nc[str_detect(tolower(all_nc), tolower(str_replace_all(pattern, "\\\\", "")))]
  # 上面这行对正则不够严格，换成真正 regex match：
  hit <- all_nc[str_detect(basename(all_nc), regex(pattern, ignore_case = TRUE))]
  if (length(hit) == 0) return(NA_character_)
  # 如果多个，优先 depthsurf
  hit2 <- hit[str_detect(basename(hit), regex("depthsurf", ignore_case = TRUE))]
  if (length(hit2) > 0) hit <- hit2
  hit[1]
}

var2file <- sapply(file_patterns, pick_file)

message("\n[FILES] matched nc files:")
print(basename(var2file))

if (any(is.na(var2file))) {
  missing <- names(var2file)[is.na(var2file)]
  stop(
    "Cannot find nc for: ", paste(missing, collapse = ", "),
    "\nUnder: ", scenario_dir,
    "\nTip: run `basename(list.files(scenario_dir, pattern='\\\\.nc$', full.names=TRUE))` and check naming."
  )
}

## -----------------------------
## 5) Extract all vars on ocean cells
## -----------------------------
# method: bilinear for continuous vars
method_map <- c(temp="bilinear", sal="bilinear", chl="bilinear", no3="bilinear", po4="bilinear", si="bilinear", o2="bilinear")

for (v in names(var2file)) {
  message("[EXTRACT] ", v, " <- ", basename(var2file[[v]]))
  grid[[v]] <- extract_on_ocean(var2file[[v]], r0, ocean_cells, method = method_map[[v]])
  gc()
}

## -----------------------------
## 6) Basic sanity checks + write
## -----------------------------
message("\n[CHECK] summaries (ocean only):")
print(sapply(grid %>% select(-Longitude, -Latitude), \(x) summary(x)[c(1,2,3,4,5,6)]))

write_tsv(grid, out_tsv)
message("\n[DONE ✅] wrote: ", normalizePath(out_tsv))
message("[ROWS] ", nrow(grid), "  [COLS] ", ncol(grid))
