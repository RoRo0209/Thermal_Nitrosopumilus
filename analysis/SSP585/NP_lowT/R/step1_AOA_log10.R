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
## AOA macroecology prep (matrix format)
## abundance_xlsx: rows=AOA genomes, cols=samples (DRR/SRR/ERR...)
## coord_xlsx: Genome=sample_id, Longitude, Latitude
## Output: AOA_sample_abundance.tsv  (ready for bioclim extraction)
## =========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
})

setwd(thermal_path("SSP585/NP_lowT"))

abundance_xlsx <- "NP_lowT_coastal_abundance.xlsx"
coord_xlsx     <- "genome_lon_lat.xlsx"
out_file       <- "NP_lowT_sample_abundance.tsv"

clean_names <- function(x) {
  x <- gsub("\ufeff", "", x, fixed = TRUE)
  x <- trimws(x)
  x
}

## 1) Read
mat <- read_excel(abundance_xlsx)
coord0 <- read_excel(coord_xlsx)

colnames(mat)   <- clean_names(colnames(mat))
colnames(coord0) <- clean_names(colnames(coord0))

stopifnot(all(c("Genome") %in% colnames(mat)))
stopifnot(all(c("Genome","Longitude","Latitude") %in% colnames(coord0)))

message("Abundance table dims: ", nrow(mat), " x ", ncol(mat))
message("Coord table dims: ", nrow(coord0), " x ", ncol(coord0))

## 2) Convert abundance matrix -> per-sample total abundance
# sample columns are all columns except the first "Genome"
sample_cols <- setdiff(colnames(mat), "Genome")

# coerce to numeric safely
mat_num <- mat
mat_num[sample_cols] <- lapply(mat_num[sample_cols], function(x) suppressWarnings(as.numeric(x)))

# sum across all AOA genomes for each sample (column-wise sum)
sample_sum <- colSums(as.data.frame(mat_num[sample_cols]), na.rm = TRUE)

aoa_sample <- tibble(
  sample_id = names(sample_sum),
  AOA_ab    = as.numeric(sample_sum)
)

# quick sanity
message("Samples with nonzero AOA_ab: ", sum(aoa_sample$AOA_ab > 0), " / ", nrow(aoa_sample))

## 3) Prepare coordinates (Genome column actually holds sample IDs)
coord <- coord0 %>%
  transmute(
    sample_id = as.character(Genome),
    Longitude = suppressWarnings(as.numeric(Longitude)),
    Latitude  = suppressWarnings(as.numeric(Latitude))
  )

## 4) Join and filter
df <- aoa_sample %>%
  left_join(coord, by = "sample_id")

if (any(is.na(df$Longitude)) || any(is.na(df$Latitude))) {
  missing <- df %>% filter(is.na(Longitude) | is.na(Latitude)) %>% pull(sample_id)
  message("WARNING: missing coordinates for (up to 20): ", paste(head(missing, 20), collapse = ", "))
}

df <- df %>%
  filter(!is.na(Longitude), !is.na(Latitude))

message("Joined points with coords: n=", nrow(df))

## 5) Handle duplicates at identical coordinates (average)
df2 <- df %>%
  group_by(Longitude, Latitude) %>%
  summarise(
    AOA_ab   = mean(AOA_ab, na.rm = TRUE),
    n_sample = n(),
    .groups = "drop"
  )

## 6) log10 transform with epsilon
positive_vals <- df2$AOA_ab[df2$AOA_ab > 0]
if (length(positive_vals) == 0) stop("All AOA_ab are <=0 after summing; check abundance values.")
eps <- min(positive_vals) / 2

df2 <- df2 %>%
  mutate(
    eps = eps,
    AOA_log10 = log10(AOA_ab + eps)
  )

write.table(df2, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", out_file)
print(head(df2, 10))
