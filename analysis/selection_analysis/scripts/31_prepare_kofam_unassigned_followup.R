#!/usr/bin/env Rscript

merged_path <- file.path(
  "selection_analysis", "results", "30_figure7c_kofam_integration",
  "Figure7c_110_families_with_Kofam.tsv"
)
fasta_path <- file.path(
  "selection_analysis", "results", "29_figure7c_kegg_submission",
  "Figure7c_110_family_representatives_BlastKOALA.faa"
)
out_dir <- file.path(
  "selection_analysis", "results", "31_kofam_unassigned_followup"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- read.delim(
  merged_path,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

read_fasta <- function(path) {
  z <- readLines(path, warn = FALSE)
  h <- which(grepl("^>", z))
  ids <- sub("^>(\\S+).*$", "\\1", z[h])
  ends <- c(h[-1] - 1L, length(z))
  seqs <- vapply(seq_along(h), function(i) {
    paste0(z[(h[i] + 1L):ends[i]], collapse = "")
  }, character(1))
  names(seqs) <- ids
  seqs
}

write_fasta <- function(ids, seqs, path, width = 70L) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  for (id in ids) {
    s <- seqs[[id]]
    if (is.null(s) || is.na(s) || !nzchar(s)) {
      stop("Missing sequence for: ", id)
    }
    writeLines(paste0(">", id), con)
    starts <- seq.int(1L, nchar(s), by = width)
    writeLines(substring(s, starts, pmin(starts + width - 1L, nchar(s))), con)
  }
}

seqs <- read_fasta(fasta_path)
if (!all(x$query_id %in% names(seqs))) {
  stop("Some merged query IDs are missing from the representative FASTA.")
}

unassigned <- x[x$Kofam_KO_count == 0L, ]
still_hypothetical <- unassigned[
  unassigned$current_annotation == "hypothetical protein",
]
selection_priority <- still_hypothetical[
  toupper(as.character(still_hypothetical$panelC_direct_support)) == "TRUE" |
    (!is.na(still_hypothetical$branchsite_foreground) &
       nzchar(trimws(still_hypothetical$branchsite_foreground))),
]

write_fasta(
  unassigned$query_id,
  seqs,
  file.path(out_dir, "Kofam_unassigned_59_all.faa")
)
write_fasta(
  still_hypothetical$query_id,
  seqs,
  file.path(out_dir, "Kofam_unassigned_42_hypothetical.faa")
)
write_fasta(
  selection_priority$query_id,
  seqs,
  file.path(out_dir, "Kofam_unassigned_13_selection_priority.faa")
)

write.table(
  unassigned,
  file.path(out_dir, "Kofam_unassigned_59_all_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  still_hypothetical,
  file.path(out_dir, "Kofam_unassigned_42_hypothetical_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  selection_priority,
  file.path(out_dir, "Kofam_unassigned_13_selection_priority_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

qc <- data.frame(
  set = c(
    "Kofam_unassigned_all",
    "Kofam_unassigned_still_hypothetical",
    "Kofam_unassigned_selection_priority"
  ),
  n_families = c(
    nrow(unassigned),
    nrow(still_hypothetical),
    nrow(selection_priority)
  )
)
write.table(
  qc,
  file.path(out_dir, "FOLLOWUP_FASTA_QC.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(qc, row.names = FALSE)
