# V4 Fig. S11. Validate rooted clades against the same-run ALE output before plotting.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5L) stop("Usage: Rscript plot_branch_events.R TREE SAME_RUN_UML_REC METADATA_XLSX VERIFIED_EVENTS_DIR OUTDIR")
if(.Platform$OS.type=="windows") invisible(Sys.setlocale("LC_CTYPE",".UTF-8"))
events_dir <- args[4]
# Bind event tables and mapping to the successful same-run extraction receipt.
if (!requireNamespace("jsonlite", quietly=TRUE)) stop("Install jsonlite before verifying the event receipt")
.receipt <- jsonlite::fromJSON(file.path(events_dir, "same_run_extraction.json"))
.hashes <- unlist(.receipt$outputs_md5)
if (!length(.hashes) || !all(c("crowns.tsv", "node_mapping.tsv", "DTLO_table_rgy.txt", "ale_numbered_species_tree.nwk") %in% names(.hashes))) stop("Incomplete extraction receipt; rerun extraction")
if (any(basename(names(.hashes)) != names(.hashes))) stop("Invalid receipt file names")
.observed <- unname(tools::md5sum(file.path(events_dir, names(.hashes))))
if (anyNA(.observed) || !identical(.observed, unname(as.character(.hashes)))) stop("Verified event files changed or are missing; rerun same-run extraction")

if(!file.exists(file.path(events_dir,"same_run_extraction.json"))) stop("Run same-run extraction first")
crowns <- read.delim(file.path(events_dir,"crowns.tsv"),colClasses="character")
stopifnot(setequal(crowns$group,c("HighT","LowT")),!anyDuplicated(crowns$group))
crowns <- crowns[match(c("HighT","LowT"),crowns$group),]
CONFIG <- list(tree_file=args[1], ale_file=args[2], metadata_file=args[3],
 dtlo_file=file.path(events_dir,"DTLO_table_rgy.txt"),
 family_files=setNames(file.path(events_dir,paste0("branch",crowns$ALE_ID,"_gene_DTLO.tsv")),crowns$ALE_ID),
 output_dir=args[5], expected_tips=177L,
 highlights=data.frame(group=crowns$group,ALE_ID=crowns$ALE_ID,colour=c("#D55E00","#0072B2")),
 show_tip_names=FALSE,max_bubble_mm=5,width_in=12,height_in=11)
.source_s <- grep("^S:",readLines(args[2],warn=FALSE),value=TRUE)
.verified_s <- paste(readLines(file.path(events_dir,"ale_numbered_species_tree.nwk"),warn=FALSE),collapse="")
if (length(.source_s)!=1L || gsub("[[:space:]]","",sub("^S:","",.source_s)) != gsub("[[:space:]]","",.verified_s)) stop("Plot ALE source differs from verified extraction")
MAKE_PLOTS <- TRUE       # FALSE: validate and export tables without plotting.

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

read_metadata <- function(path) {
  if (requireNamespace("readxl", quietly = TRUE)) {
    meta <- as.data.frame(readxl::read_excel(path, .name_repair = "minimal"))
  } else {
    need_packages("openxlsx")
    meta <- openxlsx::read.xlsx(path, check.names = FALSE)
  }
  check(all(c("gff_name", "name", "type") %in% names(meta)),
        "Metadata requires gff_name, name and type columns.")
  meta <- meta[c("gff_name", "name", "type")]
  meta[] <- lapply(meta, function(x) trimws(as.character(x)))
  ## The supplied sheet has a final category-only row without a genome ID/name.
  ## Ignore only such non-genome rows; never remove a tree tip.
  blank_id <- is.na(meta$gff_name) | !nzchar(meta$gff_name)
  blank_name <- is.na(meta$name) | !nzchar(meta$name)
  ignored <- sum(blank_id & blank_name)
  if (ignored) message("Ignoring ", ignored, " metadata row(s) without genome ID AND display name.")
  meta <- meta[!(blank_id & blank_name), , drop = FALSE]
  check_ids(meta$gff_name, "Metadata gff_name")
  check(!anyNA(meta) && all(nzchar(meta$name)) && all(nzchar(meta$type)),
        "Metadata contains missing/empty names or types. Fix these before plotting.")
  attr(meta, "ignored_non_genome_rows") <- ignored
  meta
}

validate_family_events <- function(files, map) {
  check(length(files) > 0L && all(file.exists(files)), "Missing downloaded family event table.")
  check_ids(names(files), "Downloaded family-table branch IDs")
  reference_families <- NULL
  result <- lapply(names(files), function(id) {
    d <- read.delim(files[[id]], check.names = FALSE, stringsAsFactors = FALSE)
    check(all(c("gene", "D", "T", "L", "O") %in% names(d)), "Family table requires gene/D/T/L/O columns.")
    check_ids(d$gene, paste("Family IDs in branch", id))
    if (is.null(reference_families)) reference_families <<- d$gene
    check(setequal(reference_families, d$gene), "Different gene-family sets in downloaded branch tables.")
    for (e in c("D", "T", "L", "O")) {
      check(is.numeric(d[[e]]) && all(is.finite(d[[e]]) & d[[e]] >= 0), "Invalid family event values.")
    }
    j <- match(id, map$ALE_ID)
    check(!is.na(j) && !map$isTip[j], "Downloaded family table has an unmapped internal branch.")
    sums <- colSums(d[c("D", "T", "L", "O")])
    expected <- unlist(map[j, c("D", "T", "L", "O")], use.names = FALSE)
    audit <- data.frame(ALE_ID = id, event = names(sums), n_families = nrow(d),
                        family_sum = unname(sums), aggregate_value = expected,
                        difference = unname(sums) - expected)
    check(all(abs(audit$difference) < 1e-6),
          paste("Downloaded family sums differ from DTLO aggregate at ALE", id))
    audit
  })
  do.call(rbind, result)
}

read_dtlo <- function(path) {
  raw <- read.table(path, header = TRUE, sep = "", quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
  check(!anyDuplicated(names(raw)) && "Branch" %in% names(raw), "Invalid DTLO header.")
  out <- parse_branch(raw$Branch)
  aliases <- list(D = c("D", "Duplications"), T = c("T", "Transfers"),
                  L = c("L", "Losses"), O = c("O", "Originations"))
  for (event in names(aliases)) {
    hit <- intersect(aliases[[event]], names(raw))
    check(length(hit) == 1L, paste("Missing or ambiguous DTLO column:", event))
    value <- suppressWarnings(as.numeric(raw[[hit]]))
    check(all(is.finite(value) & value >= 0), paste("Invalid/negative event values:", event))
    out[[event]] <- value
  }
  out
}

prepare_dtlo <- function(config) {
  need_packages("ape")
  paths <- unlist(config[c("tree_file", "ale_file", "dtlo_file", "metadata_file")])
  check(all(file.exists(paths)), paste("Missing input files:", paste(paths[!file.exists(paths)], collapse = ", ")))
  tree <- ape::read.tree(config$tree_file)
  ri <- clade_index(tree)
  nt <- ape::Ntip(tree)
  if (!is.null(config$expected_tips)) check(nt == config$expected_tips, "Unexpected tip count; check the input analysis.")
  ale <- read_ale_source(config$ale_file)
  check(setequal(tree$tip.label, ale$tree$tip.label), "Plotting tree and ALE S: have different genome sets.")
  check(nrow(ri) == nrow(ale$index) && setequal(ri$signature, ale$index$signature),
        "Rooted topology differs from ALE S:. Stop: do not prune, reroot, or reuse numerical node IDs.")
  map <- ri
  map$ALE_ID <- ale$index$ALE_ID[match(map$signature, ale$index$signature)]
  check_ids(map$ALE_ID, "ALE-to-R mapping")
  map$gff_name <- c(tree$tip.label, rep(NA_character_, tree$Nnode))
  expected_branch <- ale$records$Branch[match(map$ALE_ID, ale$records$ALE_ID)]

  dtlo <- read_dtlo(config$dtlo_file)
  check(nrow(dtlo) == nrow(map) && setequal(dtlo$Branch, expected_branch),
        "DTLO Branch set differs from ALE output (missing/extra branch or wrong terminal ALE ID).")
  j <- match(expected_branch, dtlo$Branch)
  check(all(map$isTip == dtlo$isTip[j]), "Terminal/internal classification mismatch.")
  map$Branch <- dtlo$Branch[j]
  for (event in c("D", "T", "L", "O")) map[[event]] <- dtlo[[event]][j]
  check(all(is.finite(as.matrix(map[c("D", "T", "L", "O")]))), "Unmapped event values.")

  meta <- read_metadata(config$metadata_file)
  k <- match(tree$tip.label, meta$gff_name)
  check(!anyNA(k), "Some tree tips have no metadata row; no tips will be silently dropped.")
  map$name <- c(meta$name[k], rep(NA_character_, tree$Nnode))
  map$type <- c(meta$type[k], rep(NA_character_, tree$Nnode))
  ## Duplicate display names are allowed because joins use unique original IDs.
  annotation <- data.frame(type = factor(meta$type[k]), row.names = tree$tip.label)
  highlights <- config$highlights
  if (nrow(highlights)) {
    h <- match(as.character(highlights$ALE_ID), map$ALE_ID)
    check(!anyNA(h) && all(!map$isTip[h]), "A highlighted ALE clade is absent or is a terminal.")
    check_ids(as.character(highlights$ALE_ID), "Highlight ALE IDs")
    highlights <- cbind(highlights, map[h, c("node", "n_descendants", "D", "T", "L", "O")])
  }
  totals <- data.frame(event = c("D", "T", "L", "O"),
    all_branches = colSums(map[c("D", "T", "L", "O")]),
    internal_branches = colSums(map[!map$isTip, c("D", "T", "L", "O")]),
    terminal_branches = colSums(map[map$isTip, c("D", "T", "L", "O")]), row.names = NULL)
  family_check <- validate_family_events(config$family_files, map)
  list(tree = tree, map = map, annotation = annotation, highlights = highlights,
       family_check = family_check,
       totals = totals, ale = ale, config = config,
       ignored_non_genome_rows = attr(meta, "ignored_non_genome_rows"))
}

export_audit <- function(x) {
  dir.create(x$config$output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(dat, name) write.csv(dat, file.path(x$config$output_dir, name),
                                           row.names = FALSE, fileEncoding = "UTF-8", na = "")
  write_csv(x$map, "ALE_to_R_node_events.csv")
  write_csv(x$highlights, "highlight_nodes_correct_events.csv")
  write_csv(x$totals, "DTLO_totals.csv")
  write_csv(x$family_check, "downloaded_family_events_verified.csv")
  ape::write.tree(x$ale$tree, file = file.path(x$config$output_dir, "ALE_species_tree_with_true_IDs.tre"))
  paths <- c(unlist(x$config[c("tree_file", "ale_file", "dtlo_file", "metadata_file")]),
             setNames(x$config$family_files, paste0("family_branch_", names(x$config$family_files))))
  write_csv(data.frame(role = names(paths), path = unname(paths),
                      md5 = unname(tools::md5sum(paths))), "input_manifest.csv")
  note <- c("# Corrected DTLO mapping", "",
    sprintf("- Tips retained: %d; internal nodes: %d; mapped event rows: %d.",
            ape::Ntip(x$tree), x$tree$Nnode, nrow(x$map)),
    "- All rooted clades match the explicit ALE S: tree one-to-one.",
    "- Terminal ALE IDs were read from the same .uml_rec branch records.",
    "- Branch coverage, duplicate IDs, metadata coverage and finite nonnegative events passed.",
    "- All input tips and internal nodes are retained, including the root record.",
    sprintf("- Ignored metadata-only rows without genome ID/name: %d; no tree tips removed.",
            x$ignored_non_genome_rows),
    "- Node means the R drawing index; ALE_ID means the ID used by the ALE output.",
    "- D/T/L/O values come from the supplied aggregated DTLO table, without rounding during mapping.",
    "- The two verified crown family tables have matching family sets and totals identical to the supplied aggregate table.",
    "- The local .uml_rec is one gene-family result, not the full aggregate event table.",
    "- These checks establish compatible numbering/topology; they cannot prove that an externally aggregated table came from the same run if all IDs are identical. Keep the original run provenance.",
    "", "## Figure", "",
    "- Four circular cladograms (branch lengths ignored), with metadata type rings.",
    "- Bubbles show internal-branch expected event counts; zero counts have no visible bubble.",
    "- Bubble area is proportional to count within each panel; scales differ across panels.",
    "- Panel subtitles report internal-branch sums only; all-branch sums are in DTLO_totals.csv.",
    "- HighT/LowT are the supplied annotations, not an inference from the event counts.",
    "- PDF/PNG rendering status is recorded separately in render_status.txt.")
  writeLines(note, file.path(x$config$output_dir, "validation_report.md"), useBytes = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(x$config$output_dir, "sessionInfo.txt"))
}

draw_four_panels <- function(x) {
  need_packages(c("ggtree", "ggplot2", "cowplot", "RColorBrewer"))
  types <- levels(x$annotation$type)
  palette <- RColorBrewer::brewer.pal(12, "Set3")
  type_colours <- setNames(if (length(types) <= 12) palette[seq_along(types)] else
                            grDevices::colorRampPalette(palette)(length(types)), types)
  titles <- c(D = "Duplications", T = "Transfers", L = "Losses", O = "Originations")
  colours <- c(D = "#E69F00", T = "#D55E00", L = "#009E73", O = "#0072B2")
  panels <- lapply(names(titles), function(event) {
    p <- ggtree::ggtree(x$tree, layout = "circular", branch.length = "none", size = 0.2)
    ## Plotting subsets are ordinary rows, not standalone tbl_tree objects.
    p$data <- as.data.frame(p$data)
    ## Explicit join using R drawing nodes; never equate node with ALE_ID.
    m <- match(p$data$node, x$map$node)
    check(!anyNA(m), "ggtree introduced nodes absent from the validated mapping.")
    p$data$value <- x$map[[event]][m]
    p$data$display_name <- x$map$name[m]
    for (i in seq_len(nrow(x$highlights))) {
      h <- x$highlights[i, ]
      p <- p + ggtree::geom_hilight(node = h$node, fill = h$colour, alpha = 0.09)
    }
    p <- p + ggtree::geom_point2(
      ggplot2::aes(subset = !isTip & value > 0, size = value),
      colour = colours[[event]], alpha = 0.8) +
      ggplot2::scale_size_area(name = "Events", max_size = x$config$max_bubble_mm) +
      ggplot2::labs(title = titles[[event]],
        subtitle = paste0("Internal total: ", format(round(sum(x$map[[event]][!x$map$isTip]), 2),
                                                     big.mark = ",", trim = TRUE))) +
      ggplot2::theme_void(base_family = "sans") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, colour = "grey35"),
        legend.title = ggplot2::element_text(size = 9), legend.text = ggplot2::element_text(size = 8),
        legend.position = "right", plot.margin = ggplot2::margin(6, 6, 6, 6))
    ## Open diamonds leave bubble sizes visible at the highlighted clade roots.
    for (i in seq_len(nrow(x$highlights))) {
      h <- x$highlights[i, ]
      position <- p$data[p$data$node == h$node, , drop = FALSE]
      p <- p + ggplot2::geom_point(data = position, ggplot2::aes(x = x, y = y),
        inherit.aes = FALSE, shape = 23, fill = NA, colour = h$colour, size = 3.2, stroke = 0.7)
    }
    p <- ggtree::gheatmap(p, x$annotation, width = 0.055, offset = 0.3,
                          colnames = FALSE, color = NA) +
      ggplot2::scale_fill_manual(values = type_colours, name = "Type", drop = FALSE,
                                labels = function(z) sub("^g__", "", z))
    if (isTRUE(x$config$show_tip_names)) {
      p <- p + ggtree::geom_tiplab(ggplot2::aes(label = display_name), size = 1.5,
                                  offset = 2, show.legend = FALSE)
    }
    ## Verify actual point-layer coordinates and area scaling, not just input joins.
    built <- ggplot2::ggplot_build(p)
    point_layer <- which(vapply(p$layers, function(z) inherits(z$geom, "GeomPointGGtree"), logical(1)))
    check(length(point_layer) == 1L, "Expected one DTLO bubble layer.")
    shown <- built$data[[point_layer]]
    expected <- p$data[!p$data$isTip & p$data$value > 0, , drop = FALSE]
    key <- function(a, b) paste(sprintf("%.10f", a), sprintf("%.10f", b), sep = ":")
    check_ids(key(expected$x, expected$y), "Internal bubble coordinates")
    jj <- match(key(expected$x, expected$y), key(shown$x, shown$y))
    check(nrow(shown) == nrow(expected) && !anyNA(jj), "Missing or misplaced event bubbles.")
    scaled <- built$plot$scales$get_scales("size")$map(expected$value)
    check(isTRUE(all.equal(shown$size[jj], scaled, check.attributes = FALSE)),
          "Bubble size differs from validated event values.")
    audit <- data.frame(event = event, R_node = expected$node,
                        ALE_ID = x$map$ALE_ID[match(expected$node, x$map$node)],
                        value = expected$value, x = expected$x, y = expected$y,
                        rendered_size = shown$size[jj])
    write.csv(audit, file.path(x$config$output_dir, paste0("plotted_bubbles_", event, ".csv")), row.names = FALSE)
    p
  })
  type_legend <- cowplot::get_legend(panels[[1]] + ggplot2::guides(size = "none") +
    ggplot2::theme(legend.position = "right", legend.key.height = grid::unit(5, "mm")))
  panels <- lapply(panels, function(p) p +
    ggplot2::guides(fill = "none", size = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme(legend.position = "bottom"))
  tree_grid <- cowplot::plot_grid(plotlist = panels, ncol = 2, labels = letters[1:4],
                                  label_size = 14, label_fontface = "bold")
  highlight_key <- cowplot::ggdraw() +
    cowplot::draw_label("Highlighted clades", x = 0.03, y = 0.9, hjust = 0, size = 10)
  for (i in seq_len(nrow(x$highlights))) {
    h <- x$highlights[i, ]
    highlight_key <- highlight_key + cowplot::draw_label(
      paste0(h$group, "  |  ALE", h$ALE_ID), x = 0.03, y = 0.9 - 0.2 * i,
      hjust = 0, size = 9, colour = h$colour)
  }
  legend_column <- cowplot::plot_grid(type_legend, highlight_key, ncol = 1, rel_heights = c(4, 1))
  figure <- cowplot::plot_grid(tree_grid, legend_column, ncol = 2, rel_widths = c(4.8, 1))
  stem <- file.path(x$config$output_dir, "DTLO_4panels_corrected")
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(paste0(stem, ".pdf"), figure, device = pdf_device,
                 width = x$config$width_in, height = x$config$height_in, bg = "white")
  ggplot2::ggsave(paste0(stem, ".png"), figure, dpi = 300,
                 width = x$config$width_in, height = x$config$height_in, bg = "white")
  if (interactive()) print(figure)
  invisible(figure)
}

run_corrected <- function(config = CONFIG, make_plots = TRUE) {
  x <- prepare_dtlo(config)
  export_audit(x)
  status_file <- file.path(config$output_dir, "render_status.txt")
  writeLines("Mapping validated. Figures not rendered in this run.", status_file)
  message(sprintf("PASS: %d tips, %d internal nodes, %d DTLO rows mapped.",
                  ape::Ntip(x$tree), x$tree$Nnode, nrow(x$map)))
  print(x$highlights, row.names = FALSE)
  if (make_plots) {
    tryCatch({
      x$figure <- draw_four_panels(x)
      writeLines("PDF and PNG exported successfully; inspect their layout before publication.", status_file)
    }, error = function(e) {
      writeLines(paste("Rendering failed:", conditionMessage(e)), status_file)
      stop(conditionMessage(e), call. = FALSE)
    })
  }
  message("Outputs: ", normalizePath(config$output_dir, winslash = "/"))
  invisible(x)
}

## Tests can load functions without running the configured analysis.
if (!isTRUE(getOption("dtlo.functions_only", FALSE))) {
  result <- run_corrected(CONFIG, make_plots = MAKE_PLOTS)
}
