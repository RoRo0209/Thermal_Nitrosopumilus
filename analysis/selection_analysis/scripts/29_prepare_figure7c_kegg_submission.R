#!/usr/bin/env Rscript

# Prepare one consensus-like representative protein per Figure 7c family for
# BlastKOALA/KAAS, using the exact protein set that passed the codeml pipeline.

project_root <- normalizePath(
  if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[1] else ".",
  winslash = "/", mustWork = TRUE
)

outdir <- file.path(
  project_root,
  "selection_analysis/results/29_figure7c_kegg_submission"
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

figure_file <- file.path(
  project_root,
  "selection_analysis/results/28_panelc_ranked_ratio_only_R/Figure7c_110_source_data.tsv"
)
branch_file <- file.path(
  project_root,
  "selection_analysis/results/16_dnds_lineage_history/source_data_branch_model.tsv"
)
branchsite_file <- file.path(
  project_root,
  "selection_analysis/results/16_dnds_lineage_history/source_data_branchsite_curated.tsv"
)
raw_dir <- file.path(
  project_root,
  "selection_analysis/results/05_codon_inputs_soft85/qc_proteins"
)
aln_dir <- file.path(
  project_root,
  paste0(
    "selection_analysis/workstation_stage_inspect/selection_analysis/results/",
    "07_alignments_gene_trees_soft85/protein_trimmed"
  )
)

read_fasta <- function(path) {
  x <- readLines(path, warn = FALSE)
  starts <- which(startsWith(x, ">"))
  if (length(starts) == 0) stop("No FASTA records in: ", path)
  ends <- c(starts[-1] - 1L, length(x))
  headers <- substring(x[starts], 2L)
  seqs <- vapply(
    seq_along(starts),
    function(i) paste0(x[(starts[i] + 1L):ends[i]], collapse = ""),
    character(1)
  )
  names(seqs) <- headers
  seqs
}

first_token <- function(x) sub(" .*$", "", x)

field_from_header <- function(header, field) {
  m <- regexec(paste0("(?:^| )", field, "=([^ ]+)"), header)
  z <- regmatches(header, m)[[1]]
  if (length(z) < 2) return("")
  z[2]
}

consensus_representative <- function(aligned, raw) {
  aln_headers <- names(aligned)
  aln_ids <- first_token(aln_headers)
  raw_headers <- names(raw)
  raw_ids <- first_token(raw_headers)
  if (anyDuplicated(aln_ids) || anyDuplicated(raw_ids)) {
    stop("Non-unique FASTA identifiers encountered")
  }

  lengths <- nchar(aligned)
  if (length(unique(lengths)) != 1) stop("Alignment sequences have unequal lengths")
  mat <- do.call(rbind, strsplit(aligned, "", fixed = TRUE))
  canonical <- strsplit("ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE)[[1]]

  consensus <- apply(mat, 2, function(col) {
    col <- col[col %in% canonical]
    if (length(col) == 0) return("-")
    tab <- sort(table(col), decreasing = TRUE)
    names(tab)[1]
  })

  scores <- vapply(seq_len(nrow(mat)), function(i) {
    comp <- mat[i, ] %in% canonical & consensus %in% canonical
    mismatch <- if (any(comp)) mean(mat[i, comp] != consensus[comp]) else 1
    gap_fraction <- mean(mat[i, ] == "-")
    ambiguous_fraction <- mean(!(mat[i, ] %in% c(canonical, "-")))
    mismatch + 0.05 * gap_fraction + 0.20 * ambiguous_fraction
  }, numeric(1))

  raw_index <- match(aln_ids, raw_ids)
  if (anyNA(raw_index)) {
    stop("Aligned identifiers missing from raw protein FASTA")
  }
  raw_lengths <- nchar(raw[raw_index])
  pick <- order(scores, -raw_lengths, aln_ids)[1]

  list(
    aligned_header = aln_headers[pick],
    sequence_id = aln_ids[pick],
    raw_header = raw_headers[raw_index[pick]],
    raw_sequence = unname(raw[raw_index[pick]]),
    consensus_distance = scores[pick]
  )
}

wrap_sequence <- function(x, width = 60L) {
  starts <- seq.int(1L, nchar(x), by = width)
  paste(substring(x, starts, pmin(starts + width - 1L, nchar(x))), collapse = "\n")
}

fig <- read.delim(figure_file, check.names = FALSE, stringsAsFactors = FALSE)
branch <- read.delim(branch_file, check.names = FALSE, stringsAsFactors = FALSE)
branchsite <- read.delim(branchsite_file, check.names = FALSE, stringsAsFactors = FALSE)

if (nrow(fig) != 110 || anyDuplicated(fig$split_family)) {
  stop("Figure source must contain 110 unique families")
}

file_map <- unique(branch[c("split_family", "file_stem")])
if (anyDuplicated(file_map$split_family)) stop("A family maps to multiple file stems")
fig <- merge(fig, file_map, by = "split_family", all.x = TRUE, sort = FALSE)
fig <- fig[order(fig$rank), ]
if (anyNA(fig$file_stem)) stop("Missing file stem for at least one family")

branch_keep <- branch[c(
  "split_family", "ecotype", "branch_FDR", "branch_direction",
  "branch_significant_FDR05"
)]
high <- branch_keep[branch_keep$ecotype == "HighT", ]
low <- branch_keep[branch_keep$ecotype == "LowT", ]
names(high)[3:5] <- c(
  "HighT_branch_FDR", "HighT_branch_direction", "HighT_branch_significant"
)
names(low)[3:5] <- c(
  "LowT_branch_FDR", "LowT_branch_direction", "LowT_branch_significant"
)
high$ecotype <- NULL
low$ecotype <- NULL
fig <- merge(fig, high, by = "split_family", all.x = TRUE, sort = FALSE)
fig <- merge(fig, low, by = "split_family", all.x = TRUE, sort = FALSE)

bs_keep <- branchsite[c(
  "split_family", "foreground", "branchsite_FDR_conservative",
  "BEB_sites_ge_0.95", "post_audit_status"
)]
names(bs_keep)[2:5] <- c(
  "branchsite_foreground", "branchsite_FDR",
  "BEB_sites_ge_0.95", "branchsite_status"
)
fig <- merge(fig, bs_keep, by = "split_family", all.x = TRUE, sort = FALSE)
fig <- fig[order(fig$rank), ]

fig$panelC_direct_support <- (
  fig$HighT_branch_significant == 1 &
    fig$LowT_branch_significant == 1 &
    fig$HighT_branch_direction == "lower" &
    fig$LowT_branch_direction == "higher"
)

records <- vector("list", nrow(fig))
fasta_lines <- character()

for (i in seq_len(nrow(fig))) {
  stem <- fig$file_stem[i]
  raw_file <- file.path(raw_dir, paste0(stem, ".faa"))
  aln_file <- file.path(aln_dir, paste0(stem, ".trim.faa"))
  if (!file.exists(raw_file)) stop("Missing raw protein file: ", raw_file)
  if (!file.exists(aln_file)) stop("Missing trimmed alignment: ", aln_file)

  raw <- read_fasta(raw_file)
  aligned <- read_fasta(aln_file)
  pick <- consensus_representative(aligned, raw)

  seq <- toupper(gsub("[*.-]", "", pick$raw_sequence))
  if (!grepl("^[ACDEFGHIKLMNPQRSTVWYBXZJUO]+$", seq)) {
    stop("Unexpected character in selected sequence for: ", fig$split_family[i])
  }

  safe_family <- gsub("[^A-Za-z0-9_]+", "_", fig$split_family[i])
  query_id <- sprintf("F7C%03d_%s", as.integer(fig$rank[i]), safe_family)
  fasta_lines <- c(fasta_lines, paste0(">", query_id), wrap_sequence(seq))

  records[[i]] <- data.frame(
    query_id = query_id,
    rank = fig$rank[i],
    split_family = fig$split_family[i],
    file_stem = stem,
    current_annotation = fig$annotation[i],
    selected_sequence_id = pick$sequence_id,
    selected_genome = field_from_header(pick$raw_header, "genome"),
    selected_locus = field_from_header(pick$raw_header, "locus"),
    selected_ecotype = field_from_header(pick$raw_header, "ecotype"),
    protein_length_aa = nchar(seq),
    consensus_distance = pick$consensus_distance,
    omega_HighT = fig$omega_HighT[i],
    omega_LowT = fig$omega_LowT[i],
    omega_ratio_HighT_LowT = fig$omega_ratio_HighT_LowT[i],
    HighT_branch_FDR = fig$HighT_branch_FDR[i],
    HighT_branch_direction = fig$HighT_branch_direction[i],
    LowT_branch_FDR = fig$LowT_branch_FDR[i],
    LowT_branch_direction = fig$LowT_branch_direction[i],
    panelC_direct_support = fig$panelC_direct_support[i],
    branchsite_foreground = ifelse(
      is.na(fig$branchsite_foreground[i]), "", fig$branchsite_foreground[i]
    ),
    branchsite_FDR = ifelse(is.na(fig$branchsite_FDR[i]), "", fig$branchsite_FDR[i]),
    BEB_sites_ge_0.95 = ifelse(
      is.na(fig$BEB_sites_ge_0.95[i]), "", fig$BEB_sites_ge_0.95[i]
    ),
    branchsite_status = ifelse(
      is.na(fig$branchsite_status[i]), "", fig$branchsite_status[i]
    ),
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, records)
if (nrow(manifest) != 110 || anyDuplicated(manifest$query_id)) {
  stop("Output manifest does not contain 110 unique query identifiers")
}

fasta_path <- file.path(outdir, "Figure7c_110_family_representatives_BlastKOALA.faa")
writeLines(fasta_lines, fasta_path, useBytes = TRUE)

write.table(
  manifest,
  file.path(outdir, "Figure7c_110_family_representative_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

direct <- manifest[manifest$panelC_direct_support, ]
write.table(
  direct,
  file.path(outdir, "Figure7c_direct_support_candidate_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

summary <- data.frame(
  metric = c(
    "families_in_figure",
    "representative_proteins_written",
    "unique_query_ids",
    "minimum_protein_length_aa",
    "maximum_protein_length_aa",
    "panelC_direct_support_genes",
    "branchsite_candidates_in_110"
  ),
  value = c(
    nrow(fig),
    nrow(manifest),
    length(unique(manifest$query_id)),
    min(manifest$protein_length_aa),
    max(manifest$protein_length_aa),
    sum(manifest$panelC_direct_support),
    sum(manifest$branchsite_foreground != "")
  )
)
write.table(
  summary,
  file.path(outdir, "PREPARATION_QC.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Wrote BlastKOALA input:", fasta_path, "\n")
cat("Families:", nrow(manifest), "\n")
cat(
  "Panel-C direct-support genes:",
  paste(direct$split_family, collapse = ", "),
  "\n"
)
