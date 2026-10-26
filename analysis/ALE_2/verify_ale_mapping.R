# Match the original IDs from the SAME ALE run by full descendant sets.
check <- function(ok, msg) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(msg, call. = FALSE)
}
need_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  check(!length(missing), paste("Missing R packages:", paste(missing, collapse = ", "),
    "\nInstall CRAN packages with install.packages(); install ggtree with BiocManager::install('ggtree')."))
}
check_ids <- function(x, what) {
  check(!anyNA(x) && all(nzchar(x)) && !anyDuplicated(x),
        paste(what, "contains missing, empty or duplicated identifiers."))
}

## Descendant-set signatures identify rooted clades, independently of R indices.
## Length-prefixed tip IDs avoid delimiter collisions; never use display names.
clade_index <- function(tree) {
  check(inherits(tree, "phylo"), "Expected one phylo tree, not multiple trees.")
  check_ids(tree$tip.label, "Tree tip labels")
  check(ape::is.rooted(tree), "Tree must be rooted; no automatic rerooting is allowed.")
  nt <- ape::Ntip(tree); nn <- tree$Nnode; total <- nt + nn
  kids <- split(tree$edge[, 2], tree$edge[, 1])
  memo <- vector("list", total)
  descendants <- function(node) {
    if (!is.null(memo[[node]])) return(memo[[node]])
    value <- if (node <= nt) tree$tip.label[node] else {
      child <- kids[[as.character(node)]]
      check(length(child) >= 2L, "Unary internal node found: clade mapping is ambiguous.")
      sort(unlist(lapply(child, descendants), use.names = FALSE), method = "radix")
    }
    memo[[node]] <<- value
    value
  }
  sets <- lapply(seq_len(total), descendants)
  sig <- vapply(sets, function(x) paste0(nchar(x, type = "bytes"), ":", x, collapse = ""), "")
  check_ids(sig, "Rooted clade signatures")
  data.frame(node = seq_len(total), isTip = seq_len(total) <= nt,
             signature = sig, n_descendants = lengths(sets),
             descendant_ids = vapply(sets, paste, "", collapse = ";"),
             stringsAsFactors = FALSE)
}

parse_branch <- function(branch) {
  branch <- trimws(as.character(branch))
  internal <- grepl("^[0-9]+$", branch)
  terminal <- grepl("^.+\\([0-9]+\\)$", branch)
  check(all(internal | terminal), paste("Unrecognized Branch:",
        paste(head(branch[!(internal | terminal)]), collapse = ", ")))
  ids <- branch
  ids[terminal] <- sub("^.*\\(([0-9]+)\\)$", "\\1", branch[terminal])
  tip <- rep(NA_character_, length(branch))
  tip[terminal] <- sub("\\([0-9]+\\)$", "", branch[terminal])
  check_ids(branch, "Branch column")
  check_ids(ids, "ALE IDs across terminal and internal branches")
  data.frame(Branch = branch, ALE_ID = ids, isTip = terminal,
             gff_name = tip, stringsAsFactors = FALSE)
}

read_ale_source <- function(path) {
  lines <- readLines(path, warn = FALSE)
  check(any(grepl("^#ALEml_undated", lines)),
        "The selected file is not an ALEml_undated .uml_rec output.")
  species <- grep("^S:[[:space:]]*", lines, value = TRUE)
  check(length(species) == 1L, "Expected exactly one S: species-tree line in .uml_rec.")
  tree <- ape::read.tree(text = sub("^S:[[:space:]]*", "", species))
  ci <- clade_index(tree)
  labs <- tree$node.label
  check(length(labs) == tree$Nnode && all(grepl("^[0-9]+$", labs)),
        "ALE S: tree must contain a numeric ID for every internal node.")
  check_ids(labs, "ALE S: internal labels")
  ## Read actual terminal ALE IDs from THIS file, not alphabetical guesses.
  records <- grep("^S_(terminal|internal)_branch[[:space:]]", lines, value = TRUE)
  fields <- strsplit(trimws(records), "[[:space:]]+")
  check(length(fields) == nrow(ci) && all(lengths(fields) >= 6L),
        "Incomplete or malformed per-branch records in .uml_rec.")
  kinds <- vapply(fields, `[`, "", 1L)
  ids <- parse_branch(vapply(fields, `[`, "", 2L))
  check(all(ids$isTip == (kinds == "S_terminal_branch")), "ALE branch kind/ID mismatch.")
  check(setequal(ids$gff_name[ids$isTip], tree$tip.label), "ALE terminal records differ from S: tips.")
  check(setequal(ids$ALE_ID[!ids$isTip], labs), "ALE internal records differ from S: labels.")
  ci$ALE_ID <- NA_character_
  nt <- ape::Ntip(tree)
  ci$ALE_ID[seq_len(nt)] <- ids$ALE_ID[match(tree$tip.label, ids$gff_name)]
  ci$ALE_ID[nt + seq_len(tree$Nnode)] <- labs
  check_ids(ci$ALE_ID, "Complete ALE mapping")
  list(tree = tree, index = ci, records = ids,
       header = lines[grep("^#ALEml_undated", lines)[1]])
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript verify_ale_mapping.R SPECIES_TREE SAME_RUN_UML_REC GROUPS_TSV OUTDIR")
need_packages("ape")
tree <- ape::read.tree(args[1]); source <- read_ale_source(args[2])
ci <- clade_index(tree)
check(setequal(ci$signature, source$index$signature), "Input and ALE rooted clades differ; do not reuse this numbering.")
ci$ALE_ID <- source$index$ALE_ID[match(ci$signature, source$index$signature)]
groups <- read.delim(args[3], colClasses = "character", check.names = FALSE)
check(all(c("genome_id", "group") %in% names(groups)), "GROUPS_TSV needs genome_id and group.")
check_ids(groups$genome_id, "Group genome IDs")
check(all(groups$genome_id %in% tree$tip.label), "Group IDs must match original tree tip IDs, not display names.")
crowns <- lapply(c("HighT", "LowT"), function(g) {
  tips <- sort(groups$genome_id[groups$group == g], method = "radix")
  check(length(tips) >= 2L, paste("Too few group members:",g))
  sig <- paste0(nchar(tips, type="bytes"), ":", tips, collapse="")
  j <- match(sig, ci$signature)
  check(!is.na(j) && !ci$isTip[j], paste("Group must match a complete rooted clade exactly:",g))
  data.frame(group=g, ALE_ID=ci$ALE_ID[j], R_node=ci$node[j], n_descendants=length(tips),
             descendant_ids=paste(tips,collapse=";"), stringsAsFactors=FALSE)
})
crowns <- do.call(rbind,crowns)
check(!anyDuplicated(crowns$ALE_ID), "The two crowns cannot be the same branch.")
dir.create(args[4], recursive=TRUE, showWarnings=FALSE)
write.table(ci,file.path(args[4],"node_mapping.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(crowns,file.path(args[4],"crowns.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
lines <- readLines(args[2],warn=FALSE)
writeLines(sub("^S:[[:space:]]*", "", grep("^S:",lines,value=TRUE)),file.path(args[4],"ale_numbered_species_tree.nwk"))
provenance <- data.frame(role=c("species_tree","same_run_uml_rec","group_definitions"),
                         file=basename(args[1:3]),md5=unname(tools::md5sum(args[1:3])))
write.table(provenance,file.path(args[4],"mapping_sources.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("PASS:",ape::Ntip(tree),"tips;",tree$Nnode,"internal nodes; all clades uniquely matched.\n")
print(crowns[c("group","ALE_ID","R_node","n_descendants")],row.names=FALSE)
