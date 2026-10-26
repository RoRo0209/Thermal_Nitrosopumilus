# Release path adapter: data and outputs live in a separate user workspace.
.release_root <- Sys.getenv("THERMAL_PROJECT_ROOT")
if (!nzchar(.release_root)) stop("Set THERMAL_PROJECT_ROOT to the analysis workspace")
.release_root <- normalizePath(.release_root, winslash = "/", mustWork = TRUE)
.release_input <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) stop("Set required input environment variable: ", name)
  normalizePath(value, winslash = "/", mustWork = TRUE)
}
setwd(.release_root)

out <- "TEMP_experiment/Figure_Temperature_Incubation"
d <- read.delim(file.path(out, "SourceData_Temperature_Incubation_summary.tsv"), check.names = FALSE)
expected <- c("Sample", "Temperature", "Time", "Replicate", "HighT_abundance", "LowT_abundance",
              "Total_two_group", "HighT_fraction", "LowT_fraction")
stopifnot(identical(names(d), expected), nrow(d) == 19, all(d[["Total_two_group"]] > 0))

fraction_error <- max(abs(d[["HighT_fraction"]] + d[["LowT_fraction"]] - 1))
total_error <- max(abs(d[["HighT_abundance"]] + d[["LowT_abundance"]] - d[["Total_two_group"]]))

h <- read.delim(file.path(out, "SourceData_Temperature_Drivers_heatmap.tsv"), check.names = FALSE)
s <- read.delim(file.path(out, "Selected_driver_genomes.tsv"), check.names = FALSE)
stopifnot(nrow(s) == 12, all(table(s[["Metadata_ecotype"]]) == 6), nrow(h) == 228)

files <- c(
  "Figure_Temperature_Incubation_2panel_equal_width.pdf",
  "Figure_Temperature_Incubation_2panel_equal_width.svg",
  "Figure_Temperature_Incubation_2panel_equal_width.png",
  "Supplementary_Fig_Temperature_Composition.pdf",
  "Supplementary_Fig_Temperature_Composition.svg",
  "Supplementary_Fig_Temperature_Composition.png",
  "Supplementary_Fig_Temperature_Drivers_heatmap.pdf",
  "Supplementary_Fig_Temperature_Drivers_heatmap.svg",
  "Supplementary_Fig_Temperature_Drivers_heatmap.png"
)
file_information <- file.info(file.path(out, files))
stopifnot(all(file_information[["size"]] > 10000))

svg_text <- paste(unlist(lapply(file.path(out, files[c(2, 5, 8)]), readLines, warn = FALSE)), collapse = "")
stopifnot(!grepl(">B0C<|>B7<", svg_text))

cat("summary rows:", nrow(d), "\n")
cat("fraction closure max error:", format(fraction_error, scientific = TRUE), "\n")
cat("abundance sum max error:", format(total_error, scientific = TRUE), "\n")
cat("heatmap rows:", nrow(h), "; selected genomes:", nrow(s), "\n")
cat("exports present and non-empty: TRUE\n")

if (requireNamespace("png", quietly = TRUE)) {
  main_png <- png::readPNG(file.path(out, files[3]), native = TRUE)
  heat_png <- png::readPNG(file.path(out, files[9]), native = TRUE)
  cat("main PNG pixels:", dim(main_png)[2], "x", dim(main_png)[1], "\n")
  cat("heatmap PNG pixels:", dim(heat_png)[2], "x", dim(heat_png)[1], "\n")
}

print(aggregate(cbind(HighT_fraction, Total_two_group) ~ Time + Temperature, d,
                function(x) c(mean = mean(x), sd = sd(x))))
