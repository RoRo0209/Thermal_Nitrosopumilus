#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript selection_analysis/scripts/32_merge_blastkoala_followup.R <BlastKOALA_user_ko.txt>")
}

blast_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
kofam_path <- file.path(
  "selection_analysis", "results", "30_figure7c_kofam_integration",
  "Figure7c_110_families_with_Kofam.tsv"
)
out_dir <- file.path(
  "selection_analysis", "results", "32_blastkoala_integration"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- read.delim(
  kofam_path,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

raw_lines <- readLines(blast_path, warn = FALSE)
raw_lines <- raw_lines[nzchar(trimws(raw_lines))]
parts <- strsplit(raw_lines, "\t", fixed = TRUE)
blast <- data.frame(
  query_id = vapply(parts, function(z) trimws(z[[1]]), character(1)),
  BlastKOALA_KO = vapply(
    parts,
    function(z) if (length(z) >= 2L) trimws(z[[2]]) else "",
    character(1)
  ),
  stringsAsFactors = FALSE
)

if (nrow(blast) != 59L || anyDuplicated(blast$query_id)) {
  stop("Expected 59 unique BlastKOALA query IDs.")
}
if (!all(blast$query_id %in% x$query_id)) {
  stop("BlastKOALA result contains IDs absent from the 110-family table.")
}
if (any(x$Kofam_KO_count[match(blast$query_id, x$query_id)] != 0L)) {
  stop("BlastKOALA follow-up should contain only Kofam-unassigned families.")
}

ko_def <- c(
  K07573 = "exosome complex component CSL4",
  K04488 = "nitrogen fixation protein NifU and related proteins (IscU/NifU family)",
  K00939 = "adenylate kinase [EC:2.7.4.3]",
  K03775 = "FKBP-type peptidyl-prolyl cis-trans isomerase SlyD [EC:5.2.1.8]",
  K02977 = "small subunit ribosomal protein S27Ae",
  K00945 = "cytidylate kinase [EC:2.7.4.25]",
  K02858 = "3,4-dihydroxy 2-butanone 4-phosphate synthase [EC:4.1.99.12]",
  K02946 = "small subunit ribosomal protein S10",
  K00030 = "isocitrate dehydrogenase (NAD+) [EC:1.1.1.41]",
  K00966 = "mannose-1-phosphate guanylyltransferase [EC:2.7.7.13]",
  K01591 = "orotidine-5'-phosphate decarboxylase [EC:4.1.1.23]",
  K00615 = "transketolase [EC:2.2.1.1]",
  K03785 = "3-dehydroquinate dehydratase I [EC:4.2.1.10]",
  K05896 = "segregation and condensation protein A",
  K00928 = "aspartate kinase [EC:2.7.2.4]",
  K02876 = "large subunit ribosomal protein L15",
  K02316 = "DNA primase DnaG",
  K02948 = "small subunit ribosomal protein S11"
)

assigned_kos <- unique(blast$BlastKOALA_KO[nzchar(blast$BlastKOALA_KO)])
if (!all(assigned_kos %in% names(ko_def))) {
  stop("Missing KO definitions for: ",
       paste(setdiff(assigned_kos, names(ko_def)), collapse = ", "))
}

blast$BlastKOALA_assigned <- nzchar(blast$BlastKOALA_KO)
blast$BlastKOALA_definition <- ifelse(
  blast$BlastKOALA_assigned,
  unname(ko_def[blast$BlastKOALA_KO]),
  ""
)

concordant_existing <- c(
  "F7C008_csl4", "F7C018_iscU2", "F7C020_adkA",
  "F7C029_group_2754", "F7C044_cmk", "F7C047_ribB",
  "F7C051_group_7924", "F7C069_pyrF", "F7C074_aroD",
  "F7C093_dnaG"
)
new_informative <- c(
  "F7C039_group_4599", "F7C079_group_3131",
  "F7C081_group_1897", "F7C092_group_2958",
  "F7C109_group_2535"
)
potential_conflict <- c(
  "F7C065_aksF", "F7C066_glmU", "F7C071_dxs_2"
)

blast$BlastKOALA_review_status <- ifelse(
  !blast$BlastKOALA_assigned,
  "unassigned",
  ifelse(
    blast$query_id %in% concordant_existing,
    "concordant_with_existing_annotation",
    ifelse(
      blast$query_id %in% new_informative,
      "new_informative_annotation",
      ifelse(
        blast$query_id %in% potential_conflict,
        "potential_paralog_or_substrate_conflict",
        "manual_review_required"
      )
    )
  )
)

review_notes <- c(
  F7C065_aksF = paste(
    "Current annotation is homoisocitrate dehydrogenase;",
    "K00030 is NAD-dependent isocitrate dehydrogenase.",
    "Retain the current family-level annotation pending phylogenetic/context review."
  ),
  F7C066_glmU = paste(
    "Current annotation is bifunctional GlmU;",
    "K00966 is mannose-1-phosphate guanylyltransferase.",
    "Do not overwrite the current annotation."
  ),
  F7C071_dxs_2 = paste(
    "Current annotation is 1-deoxy-D-xylulose-5-phosphate synthase;",
    "K00615 is transketolase.",
    "Both are thiamine-dependent enzyme families; do not overwrite without phylogeny/domain context."
  )
)
blast$BlastKOALA_review_note <- ""
idx_note <- match(blast$query_id, names(review_notes))
has_note <- !is.na(idx_note)
blast$BlastKOALA_review_note[has_note] <- unname(review_notes[idx_note[has_note]])

x <- merge(x, blast, by = "query_id", all.x = TRUE, sort = FALSE)
x <- x[match(read.delim(
  kofam_path, sep = "\t", header = TRUE, quote = "",
  comment.char = "", check.names = FALSE, stringsAsFactors = FALSE
)$query_id, x$query_id), ]

x$BlastKOALA_assigned[is.na(x$BlastKOALA_assigned)] <- FALSE
x$BlastKOALA_KO[is.na(x$BlastKOALA_KO)] <- ""
x$BlastKOALA_definition[is.na(x$BlastKOALA_definition)] <- ""
x$BlastKOALA_review_status[is.na(x$BlastKOALA_review_status)] <- "not_submitted_Kofam_assigned"
x$BlastKOALA_review_note[is.na(x$BlastKOALA_review_note)] <- ""

x$Integrated_KO <- ifelse(
  x$Kofam_KO_count > 0L,
  x$Kofam_primary_KO,
  x$BlastKOALA_KO
)
x$Integrated_KO_definition <- ifelse(
  x$Kofam_KO_count > 0L,
  x$Kofam_primary_definition,
  x$BlastKOALA_definition
)
x$Integrated_KO_source <- ifelse(
  x$Kofam_KO_count > 0L,
  "KofamKOALA",
  ifelse(x$BlastKOALA_assigned, "BlastKOALA", "none")
)
x$Integrated_KO_assigned <- nzchar(x$Integrated_KO)

conflict_flag <- x$BlastKOALA_review_status ==
  "potential_paralog_or_substrate_conflict"
x$Recommended_annotation <- ifelse(
  conflict_flag,
  x$current_annotation,
  ifelse(
    x$Integrated_KO_assigned,
    x$Integrated_KO_definition,
    x$current_annotation
  )
)
x$Recommended_annotation_evidence <- ifelse(
  x$Kofam_KO_count > 0L,
  "KofamKOALA_threshold_pass",
  ifelse(
    x$BlastKOALA_review_status == "concordant_with_existing_annotation",
    "BlastKOALA_concordant",
    ifelse(
      x$BlastKOALA_review_status == "new_informative_annotation",
      "BlastKOALA_new_informative",
      ifelse(
        conflict_flag,
        "BlastKOALA_conflict_keep_current",
        ifelse(
          x$current_annotation == "hypothetical protein",
          "hypothetical_no_KO",
          "existing_annotation_no_KO"
        )
      )
    )
  )
)

direct_flag <- toupper(as.character(x$panelC_direct_support)) == "TRUE"
branchsite_flag <- !is.na(x$branchsite_foreground) &
  nzchar(trimws(x$branchsite_foreground))

write.table(
  blast,
  file.path(out_dir, "BlastKOALA_59_followup_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  blast[blast$BlastKOALA_assigned, ],
  file.path(out_dir, "BlastKOALA_18_assignments_review.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  x,
  file.path(out_dir, "Figure7c_110_families_Kofam_BlastKOALA.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  x[direct_flag, ],
  file.path(out_dir, "Figure7c_direct_support_final_annotation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
write.table(
  x[branchsite_flag, ],
  file.path(out_dir, "Figure7c_branchsite_final_annotation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

priority_hypothetical <- (
  direct_flag | branchsite_flag
) & x$current_annotation == "hypothetical protein"

qc <- data.frame(
  metric = c(
    "figure7c_families",
    "Kofam_assigned_families",
    "BlastKOALA_followup_rows",
    "BlastKOALA_new_KO_assignments",
    "BlastKOALA_unassigned",
    "combined_families_with_KO",
    "combined_KO_coverage_percent",
    "combined_families_without_KO",
    "BlastKOALA_concordant_existing",
    "BlastKOALA_new_informative",
    "BlastKOALA_potential_conflicts",
    "direct_support_with_KO",
    "branchsite_with_KO",
    "priority_hypothetical_families",
    "priority_hypothetical_rescued_by_BlastKOALA",
    "priority_hypothetical_still_without_KO"
  ),
  value = c(
    nrow(x),
    sum(x$Kofam_KO_count > 0L),
    nrow(blast),
    sum(blast$BlastKOALA_assigned),
    sum(!blast$BlastKOALA_assigned),
    sum(x$Integrated_KO_assigned),
    round(100 * mean(x$Integrated_KO_assigned), 1),
    sum(!x$Integrated_KO_assigned),
    sum(blast$BlastKOALA_review_status == "concordant_with_existing_annotation"),
    sum(blast$BlastKOALA_review_status == "new_informative_annotation"),
    sum(blast$BlastKOALA_review_status == "potential_paralog_or_substrate_conflict"),
    sum(x$Integrated_KO_assigned[direct_flag]),
    sum(x$Integrated_KO_assigned[branchsite_flag]),
    sum(priority_hypothetical),
    sum(priority_hypothetical & x$BlastKOALA_assigned),
    sum(priority_hypothetical & !x$Integrated_KO_assigned)
  ),
  stringsAsFactors = FALSE
)

write.table(
  qc,
  file.path(out_dir, "BlastKOALA_integration_QC.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
file.copy(
  blast_path,
  file.path(out_dir, "BlastKOALA_raw_user_ko.txt"),
  overwrite = TRUE
)

cat("BlastKOALA integration complete\n")
print(qc, row.names = FALSE)
