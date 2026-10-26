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
setwd(thermal_path("1.phylogenetic"))

## =========================================================
## Supplementary Figure S1: Combined Panel A + B
## Panel A: Continuous OGT contMap
## Panel B: 4-state ASR (make.simmap, following 参考.txt)
## Data: NP_nearshore_coastal_abundance.xlsx
## =========================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(openxlsx)
  library(grDevices)
})

setwd(thermal_path("1.phylogenetic"))

## =========================================================
## 1. Load tree & metadata
## =========================================================

tree_ml <- read.tree("gtdbtk.ar53.NP.treefile")

meta <- read.xlsx("chiplot-tree.xlsx", sheet = 1)
names(meta)[names(meta) == "predOGT.(C)"] <- "OGT"
outgroup_tips <- setdiff(tree_ml$tip.label, meta$FileName)
cat("Outgroup tips:", length(outgroup_tips), "\n")

## Panel A display tree: chronos on full tree, then root on outgroup
tree_full <- chronos(tree_ml, lambda = 1, quiet = TRUE)
tree_disp <- root(tree_full, outgroup = outgroup_tips, resolve.root = TRUE)

## Panel B ASR tree: root full ML tree on outgroup first, then drop outgroup
## This gives a properly rooted NP subtree (not a trifurcation)
tree_ml_rooted <- root(tree_ml, outgroup = outgroup_tips, resolve.root = TRUE)
tree_np_ml <- drop.tip(tree_ml_rooted, outgroup_tips)
tree_asr <- chronos(tree_np_ml, lambda = 1, quiet = TRUE)
cat("ASR tree:", Ntip(tree_asr), "tips, ultrametric:", is.ultrametric(tree_asr), "\n")
cat("ASR tree rooted:", is.rooted(tree_asr), "\n")

## =========================================================
## 2. Read abundance & OGT, take 3-way intersection
##    (following 参考.txt logic)
## =========================================================

## --- Nearshore coastal relative abundance ---
coast_abund <- read.xlsx("ASR/NP_nearshore_coastal_abundance.xlsx", sheet = 1)
names(coast_abund)[1] <- "genome"

## --- OGT from metadata ---
ogt_tab <- data.frame(
  FileName = meta$FileName,
  OGT = meta$OGT,
  stringsAsFactors = FALSE
)
ogt_tab <- ogt_tab[!is.na(ogt_tab$OGT), ]

## --- 3-way intersection: tree NP tips ∩ OGT ∩ abundance ---
np_tips <- tree_asr$tip.label
common <- Reduce(intersect, list(np_tips, ogt_tab$FileName, coast_abund$genome))
cat("3-way common genomes:", length(common), "\n")

## Prune ASR tree to common set (same as 参考.txt: drop.tip non-common)
if (length(setdiff(tree_asr$tip.label, common)) > 0) {
  tree_asr <- drop.tip(tree_asr, setdiff(tree_asr$tip.label, common))
}

## Align OGT and abundance to tree order
ogt_tab <- ogt_tab[match(tree_asr$tip.label, ogt_tab$FileName), ]
coast_abund <- coast_abund[match(tree_asr$tip.label, coast_abund$genome), ]
stopifnot(all(tree_asr$tip.label == ogt_tab$FileName))
stopifnot(all(tree_asr$tip.label == coast_abund$genome))

## =========================================================
## 3. Construct Trait 1: highT / lowT (OGT median threshold)
## =========================================================

OGT_vec_asr <- setNames(ogt_tab$OGT, ogt_tab$FileName)
thr_ogt <- median(OGT_vec_asr, na.rm = TRUE)
cat("OGT threshold =", thr_ogt, "C\n")
T_bin <- ifelse(OGT_vec_asr <= thr_ogt, 0, 1)

## =========================================================
## 4. Construct Trait 2: wide / non-wide
##    (abund_thr = 1e-4, freq_thr = 0.3, same as 参考.txt)
## =========================================================

abund_df <- coast_abund[, -1]
abund_df[] <- lapply(abund_df, as.numeric)
abund_mat <- as.matrix(abund_df)
rownames(abund_mat) <- coast_abund$genome
abund_mat <- abund_mat[!is.na(rownames(abund_mat)), ]
stopifnot(storage.mode(abund_mat) == "double")

abund_thr <- 1e-4
freq_thr  <- 0.31
presence_mat <- abund_mat >= abund_thr
freq_vec <- rowSums(presence_mat) / ncol(abund_mat)
C_bin <- ifelse(freq_vec >= freq_thr, 1, 0)
names(C_bin) <- rownames(abund_mat)

cat("Wide:", sum(C_bin == 1), " Non-wide:", sum(C_bin == 0), "\n")

## =========================================================
## 5. Combine into 4-state trait (same encoding as 参考.txt)
##    1 = lowT + non-wide, 2 = lowT + wide
##    3 = highT + non-wide, 4 = highT + wide
## =========================================================

T_bin <- T_bin[tree_asr$tip.label]
C_bin <- C_bin[tree_asr$tip.label]
stopifnot(all(names(T_bin) == tree_asr$tip.label))
stopifnot(all(names(C_bin) == tree_asr$tip.label))

state4 <- rep(NA_integer_, length(T_bin))
names(state4) <- tree_asr$tip.label
state4[T_bin == 0 & C_bin == 0] <- 1
state4[T_bin == 0 & C_bin == 1] <- 2
state4[T_bin == 1 & C_bin == 0] <- 3
state4[T_bin == 1 & C_bin == 1] <- 4

cat("4-state counts:\n")
print(table(state4))

traits <- state4
stopifnot(all(tree_asr$tip.label %in% names(traits)))
traits <- traits[tree_asr$tip.label]
traits <- as.factor(traits)
levels(traits) <- as.character(1:4)

## =========================================================
## 6. ASR: make.simmap (stochastic mapping, as in 参考.txt)
## =========================================================

set.seed(123)
simmap_4state <- make.simmap(
  tree  = tree_asr,
  x     = traits,
  model = "ARD",
  nsim  = 50
)
summary_4state <- summary(simmap_4state)
node_probs_4 <- summary_4state$ace

cols4 <- setNames(c("royalblue", "skyblue", "red", "orange"), as.character(1:4))

## =========================================================
## 7. Strain name mapping
## =========================================================

name_map <- read.xlsx("../representative_NP.xlsx", sheet = 1)
name_lookup <- setNames(name_map$ncbi_strain_identifiers, name_map$FileName)
outgroup_short <- c(
  "GCF_000802205.1_ASM80220v2_genomic" = "Outgroup 1",
  "GCF_007826885.1_ASM782688v1_genomic" = "Outgroup 2",
  "GCF_001870125.1_ASM187012v1_genomic" = "Outgroup 3",
  "GCF_900696045.1_Nitrosocosmicus_franklandus_C13_genomic" = "N. franklandus C13",
  "GCF_008389435.1_ASM838943v1_genomic" = "Outgroup 5"
)
name_lookup <- c(name_lookup, outgroup_short)

rename_tips <- function(tr, lookup) {
  new_labels <- lookup[tr$tip.label]
  missing <- is.na(new_labels)
  new_labels[missing] <- tr$tip.label[missing]
  dupes <- duplicated(new_labels) | duplicated(new_labels, fromLast = TRUE)
  if (any(dupes)) {
    counts <- table(new_labels[dupes])
    for (nm in names(counts)) {
      idx <- which(new_labels == nm)
      for (j in seq_along(idx)) new_labels[idx[j]] <- paste0(nm, " (", j, ")")
    }
  }
  tr$tip.label <- new_labels
  return(tr)
}

## =========================================================
## 8. Panel A: Continuous OGT contMap (on full display tree)
## =========================================================

OGT_vec <- setNames(meta$OGT, meta$FileName)
OGT_vec <- OGT_vec[!is.na(OGT_vec)]
mean_ogt <- mean(OGT_vec, na.rm = TRUE)
for (og in outgroup_tips) OGT_vec[og] <- mean_ogt

tree_A <- rename_tips(tree_disp, name_lookup)
old_names <- names(OGT_vec)
names(OGT_vec) <- name_lookup[old_names]
na_keys <- is.na(names(OGT_vec))
if (any(na_keys)) names(OGT_vec)[na_keys] <- old_names[na_keys]

obj_OGT <- contMap(tree_A, OGT_vec, plot = FALSE, res = 200)

## =========================================================
## 9. Rename simmap tips for display (Panel B uses NP-only tree)
## =========================================================

## Rename tips in simmap objects for display
simmap_disp <- simmap_4state
for (i in seq_along(simmap_disp)) {
  simmap_disp[[i]] <- rename_tips(simmap_disp[[i]], name_lookup)
}
class(simmap_disp) <- c("multiSimmap", "multiPhylo")

## Node numbers for pie charts (same as 参考.txt)
nodes <- as.numeric(rownames(node_probs_4))

## =========================================================
## 10. Combined output PDF
## =========================================================

page_w <- 240 / 25.4
page_h <- 270 / 25.4

cairo_pdf("FigS1_ASR_combined.pdf", width = page_w, height = page_h,
          family = "sans")

layout(matrix(1:2, nrow = 1), widths = c(1, 1))
par(cex.main = 1, font.main = 2, xpd = NA)

## --- Panel A: contMap ---
par(mar = c(4.0, 0.6, 3.0, 1.2))
plot(obj_OGT,
     fsize   = 0.42,
     lwd     = 1.3,
     legend  = 0.5 * max(nodeHeights(tree_A)),
     outline = FALSE,
     type    = "phylogram",
     ftype   = "i",
     mar     = par("mar"),
     sig     = 1)
mtext("a", side = 3, line = 1.0, adj = 0, font = 2, cex = 1.4)
mtext(expression("Predicted OGT (" * degree * "C)"),
      side = 1, line = 2.6, cex = 0.75)

## --- Panel B: plotSimmap + node pies (exactly as 参考.txt) ---
par(mar = c(4.0, 0.6, 3.0, 1.2))
plotSimmap(
  simmap_disp[[1]],
  cols4,
  fsize = 0.42,
  lwd   = 1.3,
  ftype = "i"
)

nodelabels(
  node   = nodes,
  pie    = as.matrix(node_probs_4),
  piecol = cols4,
  cex    = 0.4
)

add.simmap.legend(
  colors   = cols4,
  prompt   = FALSE,
  x        = 0,
  y        = max(nodeHeights(tree_asr)) * 0.95,
  vertical = FALSE
)
mtext("b", side = 3, line = 1.0, adj = 0, font = 2, cex = 1.4)

dev.off()

cat("\n=== FigS1_ASR_combined.pdf generated ===\n")
cat("  Panel a: OGT contMap (outgroup-rooted display tree)\n")
cat("  Panel b: 4-state ASR (make.simmap, 50 sims, ARD model)\n")
cat("  Data: ASR/NP_nearshore_coastal_abundance.xlsx\n")
cat("DONE\n")
