# V4 Fig. 6b. Requires crowns and event tables extracted from the same ALE run.
# Functional-bin calculations below retain the source analysis's explicit rules.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=4L) stop("Usage: Rscript plot_functional_events.R EVENTS_DIR ANNOTATION_XLSX KO_HIERARCHY_JSON OUTDIR")
if (.Platform$OS.type=="windows") invisible(Sys.setlocale("LC_CTYPE",".UTF-8"))
BASE <- normalizePath(args[1],mustWork=TRUE)
# Bind event tables and mapping to the successful same-run extraction receipt.
if (!requireNamespace("jsonlite", quietly=TRUE)) stop("Install jsonlite before verifying the event receipt")
.receipt <- jsonlite::fromJSON(file.path(BASE, "same_run_extraction.json"))
.hashes <- unlist(.receipt$outputs_md5)
if (!length(.hashes) || !all(c("crowns.tsv", "node_mapping.tsv", "DTLO_table_rgy.txt", "ale_numbered_species_tree.nwk") %in% names(.hashes))) stop("Incomplete extraction receipt; rerun extraction")
if (any(basename(names(.hashes)) != names(.hashes))) stop("Invalid receipt file names")
.observed <- unname(tools::md5sum(file.path(BASE, names(.hashes))))
if (anyNA(.observed) || !identical(.observed, unname(as.character(.hashes)))) stop("Verified event files changed or are missing; rerun same-run extraction")

ANNOTATION <- args[2]; HIERARCHY <- args[3]; OUT <- args[4]
crowns <- read.delim(file.path(BASE,"crowns.tsv"),colClasses="character")
if (!all(c("group","ALE_ID") %in% names(crowns)) || !setequal(crowns$group,c("HighT","LowT")) || anyDuplicated(crowns$group)) stop("Invalid verified crowns table")
BRANCH <- setNames(crowns$ALE_ID,crowns$group)[c("HighT","LowT")]
FILES <- setNames(file.path(BASE,paste0("branch",BRANCH,"_gene_DTLO.tsv")),names(BRANCH))
if (!file.exists(file.path(BASE,"same_run_extraction.json"))) stop("Run extract_ale_events.py on this same run first")
EVENTS <- c("D", "T", "L", "O")
CATEGORIES <- c("Nitrogen & amino-acid", "Carbon & energy", "Transport & signaling",
                "DNA repair / replication", "Stress & protein homeostasis",
                "Annotated (other)", "Unannotated")
COLOURS <- c("#2E9D43", "#D2A000", "#D66AA5", "#4389B3", "#A8B62F", "#9B9D9E", "#E0E1E2")
names(COLOURS) <- CATEGORIES
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
check <- function(ok, msg) if (length(ok) != 1L || is.na(ok) || !ok) stop(msg, call. = FALSE)
need <- c("readxl", "jsonlite", "ggplot2", "cowplot", "svglite")
check(all(vapply(need, requireNamespace, logical(1), quietly = TRUE)), "Install required R packages first.")
write_tsv <- function(d, name) write.table(d, file.path(OUT, name), sep = "\t", quote = FALSE,
                                          row.names = FALSE, na = "", fileEncoding = "UTF-8")

## Read a row as one unique family annotation. Keep paralog suffixes such as _1.
ann <- as.data.frame(readxl::read_excel(ANNOTATION, sheet = "Sheet1", .name_repair = "minimal"))
check(all(c("genome", "KO", "KEGG_description", "annotation_source") %in% names(ann)), "Bad annotation columns.")
ann[] <- lapply(ann, function(z) trimws(as.character(z)))
check(!anyNA(ann$genome) && all(nzchar(ann$genome)) && !anyDuplicated(ann$genome), "Nonunique/missing annotation family IDs.")
check(!anyNA(ann$KO) && all(grepl("^K[0-9]{5}$", ann$KO)), "Expected one valid KO per annotation family.")

read_branch <- function(group) {
  d <- read.delim(FILES[[group]], check.names = FALSE, stringsAsFactors = FALSE)
  check(all(c("gene", EVENTS) %in% names(d)), "Bad event file columns.")
  check(!anyNA(d$gene) && !anyDuplicated(d$gene), "Duplicated/missing reconciliation family IDs.")
  check(all(grepl("^ASTRAL_species_tree_rooted[.]tre_.+[.]faa[.]clean[.]ufboot$", d$gene)),
        "Unexpected reconciliation name: do not guess family IDs.")
  d$family <- sub("[.]faa[.]clean[.]ufboot$", "", sub("^ASTRAL_species_tree_rooted[.]tre_", "", d$gene))
  check(!anyDuplicated(d$family), "Cleaning reconciliation names created duplicate family IDs.")
  check(nrow(d) == 5403L, "Expected 5403 family records per branch.")
  for (e in EVENTS) check(is.numeric(d[[e]]) && all(is.finite(d[[e]]) & d[[e]] >= 0), "Invalid event counts.")
  d$group <- group; d$ALE_branch <- BRANCH[[group]]
  d
}
by_group <- lapply(names(FILES), read_branch); names(by_group) <- names(FILES)
check(setequal(by_group$HighT$family, by_group$LowT$family), "Different family sets between crown branches.")
families <- by_group$HighT$family
j <- match(families, ann$genome)
family_ann <- data.frame(family = families, annotation_matched = !is.na(j),
  KO = ann$KO[j], KEGG_description = ann$KEGG_description[j], annotation_source = ann$annotation_source[j],
  stringsAsFactors = FALSE)
matched_kos <- unique(na.omit(family_ann$KO))

## The local full KO hierarchy covers every KO matched from the supplied workbook.
## Only parent membership is used; no extra annotation is assigned to unmatched families.
hier_tree <- jsonlite::fromJSON(HIERARCHY, simplifyVector = FALSE)
hr <- list()
walk <- function(obj, parents = character()) {
  label <- obj$name
  if (grepl("^K[0-9]{5}([[:space:]]|$)", label)) {
    ko <- substr(label, 1, 6)
    if (ko %in% matched_kos) hr[[length(hr) + 1L]] <<- data.frame(
      KO = ko, parent_path = paste(parents, collapse = " > "), KO_name = label,
      stringsAsFactors = FALSE)
  }
  if (length(obj$children)) for (child in obj$children) walk(child, c(parents, label))
}
walk(hier_tree)
hier <- unique(do.call(rbind, hr))
check(all(matched_kos %in% hier$KO), "Some annotation KOs are missing from the reference hierarchy.")
write_tsv(hier, "matched_KO_hierarchy_paths.tsv")

## Explicit exclusive assignment. Specific molecular functions take precedence
## over broad metabolic pathway membership, which can be nonexclusive in KEGG.
## Full candidate memberships and conflict resolution are exported for review.
rules <- data.frame(priority = seq_len(6), code = c("DNA", "Stress", "Energy_ATP", "Transport", "Nitrogen", "Carbon"),
  category = CATEGORIES[c(4, 5, 2, 3, 1, 2)],
  evidence = c(
    "KEGG 09124 replication/repair or BR03032/BR03400",
    "BR03110 chaperones; BR/PATH03050/03051 proteasome; explicit stress/proteostasis description",
    "Description identifies ATP synthase, or V/A/F-type H+/Na+-transporting ATPase machinery",
    "KEGG09131 membrane transport; BR02000 transporters; direct sensor/regulator signaling description (not downstream membership in a signaling pathway)",
    "KEGG09105/09106 amino-acid metabolism or PATH00910 nitrogen metabolism",
    "KEGG09101 carbohydrate, 09102 energy (except nitrogen), or 09103 lipid metabolism"),
  stringsAsFactors = FALSE)
write_tsv(rules, "functional_bin_rules.tsv")

stress_pattern <- paste(c("chaperon", "heat[- ]shock", "proteasom", "peroxiredoxin", "peroxidase", "catalase",
  "superoxide dismutase", "thioredoxin", "glutaredoxin", "protein disulfide", "methionine.sulfoxide reductase",
  "osmoprotect", "ectoine", "ATP.dependent.*(Clp|Lon|FtsH|Hsl)", "(Clp|Lon|FtsH|Hsl).*protease",
  "protein repair"), collapse = "|")
transport_pattern <- "sensor histidine kinase|response regulator|two.component|chemotaxis|cyclic.di.GMP|cyclic.di.AMP"
amino_primary <- paste(c("amino.?acid", "aminotransferase", "glutam", "glutamin", "asparag", "aspart",
  "arginin", "ornithin", "lysine", "lysyl", "leucine", "isoleucine", "valine", "histidin", "histidinol",
  "tryptophan", "tyrosine", "phenylalanine", "cystein", "cystathionine", "methionine", "homocystein",
  "homoserine", "serine", "threonine", "glycine", "alanine", "proline", "isopropylmalate", "acetolactate",
  "acetohydroxy", "anthranilate", "chorismate", "pyrroline.5.carboxylate", "imidazoleglycerol.phosphate"), collapse = "|")
carbon_primary <- paste(c("malate dehydrogenase", "glycerate", "glycerol", "glyceraldehyde", "phosphoglycerate",
  "glucose", "glucokinase", "fructose", "pyruvate", "lactate dehydrogenase", "citrate", "aconitate", "aconitase",
  "succinate", "succinyl.CoA", "fumarase", "fumarate", "alcohol dehydrogenase", "acetyl.CoA C.acetyltransferase",
  "sulfate adenylyltransferase", "sulfurtransferase"), collapse = "|")
ko_groups <- split(hier, hier$KO)
classify <- function(ko, description) {
  if (is.na(ko)) return(c(category = "Unannotated", reason = "No exact family match in supplied annotation workbook",
                         candidates = "Unannotated", selected_rule = "none", overlap = "FALSE"))
  hp <- paste(ko_groups[[ko]]$parent_path, collapse = " | ")
  desc <- if (is.na(description) || !nzchar(description)) ko_groups[[ko]]$KO_name[1] else description
  test <- c(
    DNA = grepl("09124 |\\[BR:ko03032\\]|\\[BR:ko03400\\]", hp),
    Stress = grepl("\\[BR:ko03110\\]|\\[(BR|PATH):ko0305[01]\\]", hp) || grepl(stress_pattern, desc, ignore.case = TRUE),
    Energy_ATP = grepl("ATP synthase|[VAF](/A)?-type H\\+/Na\\+-transporting ATPase|F-type H\\+-transporting ATPase", desc, ignore.case = TRUE),
    Transport = grepl("09131 |\\[BR:ko02000\\]", hp) || grepl(transport_pattern, desc, ignore.case = TRUE),
    Nitrogen = grepl("0910[56] |\\[PATH:ko00910\\]", hp),
    Carbon = grepl("09101 |09103 |09102 ", hp) && !grepl("\\[PATH:ko00910\\]", hp)
  )
  ix <- which(test)
  if (!length(ix)) return(c(category = "Annotated (other)", reason = "KO annotation outside the five named functional bins",
                             candidates = "Annotated (other)", selected_rule = "other", overlap = "FALSE"))
  cats <- unique(rules$category[ix])
  ## Shared metabolic pathway membership alone does not identify a primary role.
  ## Resolve only explicit molecular descriptions; retain unclear cases in 'other'.
  if (all(names(test)[ix] %in% c("Nitrogen", "Carbon")) && length(cats) > 1L) {
    amino <- grepl(amino_primary, desc, ignore.case = TRUE)
    carbon <- grepl(carbon_primary, desc, ignore.case = TRUE)
    chosen <- if (amino && !carbon) CATEGORIES[1] else if (carbon && !amino) CATEGORIES[2] else CATEGORIES[6]
    return(c(category = chosen, reason = if (chosen == CATEGORIES[6])
      "Shared metabolic membership without an unambiguous primary function: retain in annotated-other" else
      "Shared metabolic membership resolved by explicit amino-acid or central carbon/energy molecular description",
      candidates = paste(cats, collapse = " | "), selected_rule = "metabolic_overlap_review", overlap = "TRUE"))
  }
  c(category = rules$category[ix[1]], reason = rules$evidence[ix[1]],
    candidates = paste(cats, collapse = " | "), selected_rule = rules$code[ix[1]], overlap = as.character(length(cats) > 1))
}
assigned <- t(mapply(classify, family_ann$KO, family_ann$KEGG_description, USE.NAMES = FALSE))
family_ann <- cbind(family_ann, as.data.frame(assigned, stringsAsFactors = FALSE))
check(nrow(family_ann) == 5403L && !anyDuplicated(family_ann$family), "Family annotation join expanded rows.")
check(all(family_ann$category %in% CATEGORIES), "Unclassified family.")
write_tsv(family_ann, "family_functional_assignments.tsv")
write_tsv(family_ann[family_ann$overlap == "TRUE", ], "families_with_multiple_functional_candidates.tsv")
write_tsv(family_ann[!family_ann$annotation_matched, ], "families_without_workbook_annotation.tsv")
write_tsv(as.data.frame(table(factor(family_ann$category, levels = CATEGORIES))), "family_counts_by_category.tsv")

data <- do.call(rbind, by_group)
k <- match(data$family, family_ann$family)
data <- cbind(data, family_ann[k, setdiff(names(family_ann), "family"), drop = FALSE])
rownames(data) <- NULL
check(nrow(data) == 10806L && !anyDuplicated(paste(data$group, data$family)), "Event join duplicated family records.")
write_tsv(data, "source_data_family_events_with_annotations.tsv")
long <- do.call(rbind, lapply(EVENTS, function(e) data.frame(
  group = data$group, ALE_branch = data$ALE_branch, family = data$family,
  category = data$category, event = e, value = data[[e]], stringsAsFactors = FALSE)))
summary <- aggregate(value ~ group + ALE_branch + category + event, long, sum)
grid <- expand.grid(group = names(FILES), category = CATEGORIES, event = EVENTS, stringsAsFactors = FALSE)
grid$ALE_branch <- unname(BRANCH[grid$group])
summary <- merge(grid, summary, by = c("group", "ALE_branch", "category", "event"), all.x = TRUE, sort = FALSE)
summary$value[is.na(summary$value)] <- 0
totals <- aggregate(value ~ group + ALE_branch + event, long, sum)
reference <- read.delim(file.path(BASE, "DTLO_table_rgy.txt"), check.names = FALSE, stringsAsFactors = FALSE)
checks <- do.call(rbind, lapply(seq_len(nrow(totals)), function(i) {
  row <- totals[i, ]
  reference_value <- reference[as.character(reference$Branch) == row$ALE_branch, row$event]
  check(length(reference_value) == 1, "Missing or duplicated crown branch in aggregate table.")
  stack_value <- sum(summary$value[summary$group == row$group & summary$event == row$event])
  data.frame(group = row$group, ALE_branch = row$ALE_branch, event = row$event,
             all_family_sum = row$value, functional_stack_sum = stack_value, reference_value = reference_value,
             pass = abs(row$value - reference_value) < 1e-6 && abs(stack_value - reference_value) < 1e-6)
}))
check(all(checks$pass), "Family sums, functional stacks and branch aggregate values differ.")
write_tsv(checks, "event_total_validation.tsv")

ratios <- do.call(rbind, lapply(EVENTS, function(e) {
  h <- totals$value[totals$group == "HighT" & totals$event == e]
  l <- totals$value[totals$group == "LowT" & totals$event == e]
  check(l > 0, "Cannot form HighT/LowT ratio with a zero denominator.")
  data.frame(event = e, HighT = h, LowT = l, HighT_div_LowT = h / l)
}))
write_tsv(summary, "panelB_functional_event_summary.tsv")
write_tsv(ratios, "panelB_totals_and_ratios.tsv")
print(ratios, row.names = FALSE)
cat("Annotated families:", sum(family_ann$annotation_matched), "/ 5403\n")
cat("Multiple candidate bins:", sum(family_ann$overlap == "TRUE"), "\n")
print(table(factor(family_ann$category, levels = CATEGORIES)))

## Quantitative grid: four independent, zero-based y axes, two stacked bars each.
event_labels <- c(D = "Duplications", T = "Transfers", L = "Losses", O = "Originations")
summary$event <- factor(summary$event, levels = EVENTS, labels = unname(event_labels))
totals$event <- factor(totals$event, levels = EVENTS, labels = unname(event_labels))
ratios$event <- factor(ratios$event, levels = EVENTS, labels = unname(event_labels))
summary$group <- factor(summary$group, levels = c("HighT", "LowT"))
totals$group <- factor(totals$group, levels = c("HighT", "LowT"))
summary$category <- factor(summary$category, levels = CATEGORIES)
ratios$label <- sprintf("H/L = %.2f\u00d7", ratios$HighT_div_LowT)
ratios$direction_colour <- ifelse(ratios$HighT_div_LowT >= 1, "#B74E31", "#356D93")
legend_labels <- c("Nitrogen &\namino-acid", "Carbon & energy", "Transport &\nsignaling",
                   "DNA repair /\nreplication", "Stress & protein\nhomeostasis", "Annotated (other)", "Unannotated")
names(legend_labels) <- CATEGORIES

p <- ggplot2::ggplot(summary, ggplot2::aes(x = group, y = value, fill = category)) +
  ggplot2::geom_col(width = 0.64, position = ggplot2::position_stack(reverse = TRUE)) +
  ggplot2::geom_text(data = totals, ggplot2::aes(x = group, y = value, label = sprintf("%.2f", value)),
                    inherit.aes = FALSE, vjust = -0.4, size = 3.0, fontface = "bold", colour = "#333333") +
  ggplot2::geom_text(data = ratios, ggplot2::aes(x = Inf, y = Inf, label = label, colour = direction_colour),
                    inherit.aes = FALSE, hjust = 1.15, vjust = 1.45, size = 2.8, fontface = "bold") +
  ggplot2::facet_wrap(~event, ncol = 2, scales = "free_y", axes = "all_x") +
  ggplot2::scale_fill_manual(values = COLOURS, labels = legend_labels, drop = FALSE, name = NULL) +
  ggplot2::scale_colour_identity(guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.19)),
                             breaks = scales::breaks_pretty(n = 4)) +
  ggplot2::labs(x = NULL, y = "Expected event count",
                caption = paste0("HighT: ALE",BRANCH[["HighT"]],"  |  LowT: ALE",BRANCH[["LowT"]],"  |  H/L: HighT / LowT")) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3, byrow = TRUE)) +
  ggplot2::theme_classic(base_size = 10, base_family = "sans") +
  ggplot2::theme(axis.line = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA, colour = "#555555", linewidth = 0.35),
    axis.ticks = ggplot2::element_line(colour = "#555555", linewidth = 0.3),
    axis.text = ggplot2::element_text(colour = "#333333", size = 8.5),
    axis.title.y = ggplot2::element_text(size = 10, margin = ggplot2::margin(r = 7)),
    strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 11),
    panel.spacing = grid::unit(6, "mm"), legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 7.7), legend.key.size = grid::unit(3.2, "mm"),
    legend.spacing.y = grid::unit(1.4, "mm"), plot.caption = ggplot2::element_text(size = 7.6, hjust = 0),
    plot.margin = ggplot2::margin(9, 8, 4, 7))

## Verify actual rendered stack heights and one bar per group per panel.
built <- ggplot2::ggplot_build(p)
bars <- built$data[[1]]
drawn <- aggregate(ymax ~ PANEL + x, bars, max)
check(nrow(drawn) == 8L, "Expected eight total stacked bars.")
layout <- built$layout$layout
for (i in seq_len(nrow(drawn))) {
  label <- as.character(layout$event[match(as.integer(drawn$PANEL[i]), layout$PANEL)])
  group <- c("HighT", "LowT")[as.integer(drawn$x[i])]
  expected <- totals$value[as.character(totals$event) == label & totals$group == group]
  check(length(expected) == 1 && abs(drawn$ymax[i] - expected) < 1e-6, "Rendered bar height is incorrect.")
}
figure <- cowplot::ggdraw(p) + cowplot::draw_label("b", x = 0.006, y = 0.994,
                                                hjust = 0, vjust = 1, size = 14, fontface = "bold")
stem <- file.path(OUT, "PanelB_functional_events")
ggplot2::ggsave(paste0(stem, ".pdf"), figure, device = grDevices::cairo_pdf,
               width = 160, height = 155, units = "mm", bg = "white")
ggplot2::ggsave(paste0(stem, ".svg"), figure, device = svglite::svglite,
               width = 160, height = 155, units = "mm", bg = "white")
ggplot2::ggsave(paste0(stem, ".png"), figure, width = 160, height = 155, units = "mm", dpi = 600, bg = "white")
ggplot2::ggsave(file.path(OUT, "PanelB_preview.png"), figure, width = 160, height = 155, units = "mm", dpi = 180, bg = "white")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
input_files <- c(FILES, annotation = ANNOTATION, KEGG_hierarchy = HIERARCHY,
                 aggregate_table = file.path(BASE, "DTLO_table_rgy.txt"))
write_tsv(data.frame(role = names(input_files), path = unname(input_files),
                     md5 = unname(tools::md5sum(input_files))), "input_manifest.tsv")
writeLines(c("PASS: 5403 unique families per group, same family set.",
             "PASS: exact annotation join; 10806 group-family records retained.",
             "PASS: every family assigned to exactly one display category.",
             "PASS: all 8 family sums and category-stack sums match the correct ALE aggregate rows.",
             "PASS: actual rendered heights of all 8 bars match expected counts.",
             "PDF, SVG and PNG exports completed."), file.path(OUT, "validation_checks.txt"))
cat("PANEL B EXPORT AND ALL CHECKS COMPLETE\n")
