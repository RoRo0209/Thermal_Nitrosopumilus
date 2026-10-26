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
setwd(thermal_path("SSP585/NP_highT"))

## =========================================================
## AOA macroecology V0 (Bio-ORACLE baseline) analysis script
## - Input : AOA_points_with_BioORACLE.tsv (points + predictors)
## - Output: model table, VIF table, RF model, importance, PDP plots
## =========================================================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

## -----------------------------
## 0) Paths (EDIT)
## -----------------------------
setwd(thermal_path("SSP585/NP_highT"))

in_file <- "AOA_points_with_BioORACLE_nearshore_filtered.tsv"  # CHANGED: use nearshore-filtered data
out_dir <- "AOA_RF_V0_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------
## 1) Read + sanity checks
## -----------------------------
dat0 <- read_tsv(in_file, show_col_types = FALSE)
stopifnot(all(c("Longitude","Latitude") %in% names(dat0)))

# response column
y_col <- if ("AOA_log10" %in% names(dat0)) "AOA_log10" else if ("AOA_ab" %in% names(dat0)) "AOA_ab" else NA
if (is.na(y_col)) stop("Cannot find response column AOA_log10 or AOA_ab in input file.")

## -----------------------------
## 2) Detect Bio-ORACLE columns (AUTO)
##    Unify to temp/sal/o2/no3/po4/si/chl
## -----------------------------
# two possible naming styles:
# A) short names: temp sal o2 no3 po4 si chl
# B) dataset_id names: thetao_baseline... / so_baseline... etc

map_candidates <- list(
  temp = c(
    "temp",
    "thetao_baseline_2000_2019_depthsurf",
    "thetao_baseline_2000_2018_depthsurf"
  ),
  sal  = c(
    "sal",
    "so_baseline_2000_2019_depthsurf",
    "so_baseline_2000_2018_depthsurf"
  ),
  o2   = c(
    "o2",
    "o2_baseline_2000_2018_depthsurf"
  ),
  no3  = c(
    "no3",
    "no3_baseline_2000_2018_depthsurf"
  ),
  po4  = c(
    "po4",
    "po4_baseline_2000_2018_depthsurf"
  ),
  si   = c(
    "si",
    "si_baseline_2000_2018_depthsurf"
  ),
  chl  = c(
    "chl",
    "chl_baseline_2000_2018_depthsurf"
  )
)


# find which columns exist
col_map <- sapply(names(map_candidates), function(v) {
  hits <- intersect(map_candidates[[v]], names(dat0))
  if (length(hits) == 0) NA_character_ else hits[1]
})

message("Detected predictor columns:")
print(col_map)

# must have temp + sal at minimum
if (is.na(col_map["temp"]) || is.na(col_map["sal"])) {
  stop("Missing required predictors temp and/or sal.\n",
       "Current columns in file include (first 40):\n",
       paste(head(names(dat0), 40), collapse = ", "),
       "\n\nTip: run `colnames(dat0)` and check whether your extraction script produced temp/sal or thetao/so dataset_id columns.")
}

# build modeling table with unified names
dat <- dat0 %>%
  transmute(
    Longitude, Latitude,
    y   = as.numeric(.data[[y_col]]),
    temp = as.numeric(.data[[col_map["temp"]]]),
    sal  = as.numeric(.data[[col_map["sal"]]]),
    o2   = if (!is.na(col_map["o2"]))  as.numeric(.data[[col_map["o2"]]])  else NA_real_,
    no3  = if (!is.na(col_map["no3"])) as.numeric(.data[[col_map["no3"]]]) else NA_real_,
    po4  = if (!is.na(col_map["po4"])) as.numeric(.data[[col_map["po4"]]]) else NA_real_,
    si   = if (!is.na(col_map["si"]))  as.numeric(.data[[col_map["si"]]])  else NA_real_,
    chl  = if (!is.na(col_map["chl"])) as.numeric(.data[[col_map["chl"]]]) else NA_real_
  )

# quick missingness report
message("\nNA fraction by column:")
print(round(sapply(dat, function(x) mean(is.na(x))), 3))

## -----------------------------
## 3) Filter (minimal): require y + temp + sal
## -----------------------------
dat2 <- dat %>% filter(!is.na(y), !is.na(temp), !is.na(sal))
message("\nRows after filtering (need y+temp+sal): ", nrow(dat2))
write_tsv(dat2, file.path(out_dir, "AOA_model_table_V0.tsv"))

if (nrow(dat2) < 30) {
  message("\nWARNING: very few rows left (", nrow(dat2), ").",
          "\nThis is NOT a modeling issue; it means your Bio-ORACLE extraction still has lots of NA.",
          "\nFix extraction first (most points might be on land cells / coastline NA / bbox mismatch).")
}

## -----------------------------
## 4) (optional) build a smaller predictor set with enough non-NA
## -----------------------------
pred_cols <- c("temp","sal","o2","no3","po4","si","chl")
keep_pred <- pred_cols[sapply(dat2[pred_cols], function(x) sum(!is.na(x)) >= 20)]  # require >=20 non-NA
message("\nPredictors kept (>=20 non-NA): ", paste(keep_pred, collapse=", "))

# save final modeling frame
model_df <- dat2 %>% select(y, all_of(keep_pred), Longitude, Latitude)
write_tsv(model_df, file.path(out_dir, "AOA_model_table_V0_final.tsv"))

message("\nDone. Next: fit RF / GAM / PGLS using AOA_model_table_V0_final.tsv")
