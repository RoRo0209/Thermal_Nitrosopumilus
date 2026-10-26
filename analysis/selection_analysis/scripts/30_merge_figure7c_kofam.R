#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript selection_analysis/scripts/30_merge_figure7c_kofam.R <Kofam_result.txt>")
}

raw_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
manifest_path <- file.path(
  "selection_analysis", "results", "29_figure7c_kegg_submission",
  "Figure7c_110_family_representative_manifest.tsv"
)
out_dir <- file.path(
  "selection_analysis", "results", "30_figure7c_kofam_integration"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(
  manifest_path,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(manifest) != 110L || anyDuplicated(manifest$query_id)) {
  stop("The Figure 7c manifest must contain 110 unique query IDs.")
}

raw_lines <- readLines(raw_path, warn = FALSE)
data_lines <- raw_lines[
  nzchar(trimws(raw_lines)) &
    !grepl("^\\s*#", raw_lines)
]

pattern <- paste0(
  "^\\s*(\\*?)\\s*(\\S+)\\s+(K[0-9]{5})\\s+",
  "([0-9.]+)\\s+([0-9.]+)\\s+(\\S+)\\s+(.+?)\\s*$"
)

parsed <- lapply(data_lines, function(x) {
  hit <- regmatches(x, regexec(pattern, x, perl = TRUE))[[1]]
  if (length(hit) != 8L) {
    stop("Could not parse KofamKOALA line: ", x)
  }
  data.frame(
    threshold_pass = identical(hit[[2]], "*"),
    query_id = hit[[3]],
    KO = hit[[4]],
    threshold = as.numeric(hit[[5]]),
    score = as.numeric(hit[[6]]),
    evalue = hit[[7]],
    definition = hit[[8]],
    stringsAsFactors = FALSE
  )
})

hits <- do.call(rbind, parsed)
hits$score_threshold_ratio <- hits$score / hits$threshold
hits <- hits[order(hits$query_id, -hits$score_threshold_ratio), ]
row.names(hits) <- NULL

unknown_ids <- setdiff(unique(hits$query_id), manifest$query_id)
if (length(unknown_ids)) {
  stop("Kofam result contains IDs absent from the manifest: ",
       paste(unknown_ids, collapse = ", "))
}

collapse_field <- function(x) paste(x, collapse = "; ")

hit_split <- split(hits, hits$query_id)
summary_rows <- lapply(manifest$query_id, function(id) {
  z <- hit_split[[id]]
  if (is.null(z)) {
    return(data.frame(
      query_id = id,
      Kofam_status = "unassigned",
      Kofam_KO_count = 0L,
      Kofam_primary_KO = NA_character_,
      Kofam_primary_definition = NA_character_,
      Kofam_primary_score_threshold_ratio = NA_real_,
      Kofam_KOs = NA_character_,
      Kofam_definitions = NA_character_,
      Kofam_scores = NA_character_,
      Kofam_thresholds = NA_character_,
      Kofam_Evalues = NA_character_,
      Kofam_all_threshold_pass = NA,
      stringsAsFactors = FALSE
    ))
  }

  z <- z[order(-z$score_threshold_ratio), ]
  data.frame(
    query_id = id,
    Kofam_status = if (nrow(z) == 1L) "single_KO" else "multiple_KOs",
    Kofam_KO_count = nrow(z),
    Kofam_primary_KO = z$KO[[1]],
    Kofam_primary_definition = z$definition[[1]],
    Kofam_primary_score_threshold_ratio = z$score_threshold_ratio[[1]],
    Kofam_KOs = collapse_field(z$KO),
    Kofam_definitions = collapse_field(z$definition),
    Kofam_scores = collapse_field(format(z$score, trim = TRUE, scientific = FALSE)),
    Kofam_thresholds = collapse_field(format(z$threshold, trim = TRUE, scientific = FALSE)),
    Kofam_Evalues = collapse_field(z$evalue),
    Kofam_all_threshold_pass = all(z$threshold_pass),
    stringsAsFactors = FALSE
  )
})

kofam_summary <- do.call(rbind, summary_rows)
merged <- merge(
  manifest,
  kofam_summary,
  by = "query_id",
  all.x = TRUE,
  sort = FALSE
)
merged <- merged[match(manifest$query_id, merged$query_id), ]

direct_flag <- toupper(as.character(merged$panelC_direct_support)) == "TRUE"
branchsite_flag <- !is.na(merged$branchsite_foreground) &
  nzchar(trimws(merged$branchsite_foreground))

direct <- merged[direct_flag, ]
branchsite <- merged[branchsite_flag, ]

write.table(
  hits,
  file.path(out_dir, "Kofam_assignments_long.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  merged,
  file.path(out_dir, "Figure7c_110_families_with_Kofam.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  direct,
  file.path(out_dir, "Figure7c_direct_support_candidates_with_Kofam.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  branchsite,
  file.path(out_dir, "Figure7c_branchsite_candidates_with_Kofam.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

qc <- data.frame(
  metric = c(
    "figure7c_families",
    "kofam_assignment_rows",
    "families_with_at_least_one_KO",
    "families_without_KO",
    "family_annotation_rate_percent",
    "families_with_multiple_KOs",
    "all_assignment_rows_above_KO_threshold",
    "direct_support_candidates",
    "direct_support_candidates_with_KO",
    "branchsite_candidates",
    "branchsite_candidates_with_KO"
  ),
  value = c(
    nrow(merged),
    nrow(hits),
    sum(merged$Kofam_KO_count > 0L),
    sum(merged$Kofam_KO_count == 0L),
    round(100 * mean(merged$Kofam_KO_count > 0L), 1),
    sum(merged$Kofam_KO_count > 1L),
    all(hits$threshold_pass),
    nrow(direct),
    sum(direct$Kofam_KO_count > 0L),
    nrow(branchsite),
    sum(branchsite$Kofam_KO_count > 0L)
  ),
  stringsAsFactors = FALSE
)

write.table(
  qc,
  file.path(out_dir, "Kofam_integration_QC.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

file.copy(
  raw_path,
  file.path(out_dir, "KofamKOALA_raw_result.txt"),
  overwrite = TRUE
)

cat("KofamKOALA integration complete\n")
print(qc, row.names = FALSE)
