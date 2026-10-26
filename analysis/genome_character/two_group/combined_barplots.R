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
setwd(thermal_path("genome_character/two_group"))

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# ── Shared: group info ─────────────────────────────────────────────────────────
meta <- read_excel("two_group.xlsx")[, c("FileName", "group")]
meta$group <- factor(meta$group, levels = c("lowT", "highT"))
fill_scale <- scale_fill_manual(values = c(lowT = "#8EC6E8", highT = "#F4A261"),
                                name = "Group")

sig_from_padj <- function(p) ifelse(is.na(p), "ns",
                             ifelse(p < 0.001, "***",
                             ifelse(p < 0.01,  "**",
                             ifelse(p < 0.05,  "*", "ns"))))

make_wilcox_stats <- function(counts, key) {
  counts |>
    group_by(.data[[key]]) |>
    summarise(
      mean_lowT  = mean(count[group == "lowT"]),
      mean_highT = mean(count[group == "highT"]),
      p_value    = tryCatch(
        wilcox.test(count[group == "highT"], count[group == "lowT"],
                    exact = FALSE)$p.value,
        error = function(e) NA_real_),
      .groups = "drop"
    ) |>
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           sig   = sig_from_padj(p_adj))
}

# ── CAZyme ─────────────────────────────────────────────────────────────────────
caz <- read.table("overview.txt", sep = "\t", header = TRUE,
                  quote = "", comment.char = "", stringsAsFactors = FALSE,
                  check.names = FALSE)
colnames(caz) <- trimws(colnames(caz))
caz <- caz[caz[["#ofTools"]] >= 2, ]
caz$genome <- sub("__.*$", "", caz[["Gene ID"]])
caz$family <- str_extract(caz[["HMMER"]], "^(CBM|GH|GT|CE|AA|PL)")
caz <- caz[!is.na(caz$family), ]

caz_counts <- caz |> group_by(genome, family) |>
  summarise(count = n(), .groups = "drop")
caz_grid <- expand.grid(genome = unique(meta$FileName),
                        family = unique(caz_counts$family),
                        stringsAsFactors = FALSE)
caz_counts <- left_join(caz_grid, caz_counts, by = c("genome", "family"))
caz_counts$count[is.na(caz_counts$count)] <- 0
caz_counts <- left_join(caz_counts, meta, by = c("genome" = "FileName"))
caz_counts <- caz_counts[!is.na(caz_counts$group), ]
caz_stats <- make_wilcox_stats(caz_counts, "family")

caz_plot_df <- caz_counts |> group_by(family, group) |>
  summarise(mean = mean(count), se = sd(count)/sqrt(n()), .groups = "drop")
caz_sig <- left_join(
  caz_plot_df |> filter(group == "highT") |> select(family, mean),
  caz_stats |> select(family, sig), by = "family")

p_caz <- ggplot(caz_plot_df, aes(x = family, y = mean, fill = group)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_errorbar(aes(ymin = pmax(mean - se, 0), ymax = mean + se),
                position = position_dodge(0.7), width = 0.25, linewidth = 0.5) +
  geom_text(data = caz_sig, aes(x = family, y = mean + 0.5, label = sig),
            inherit.aes = FALSE, size = 4, vjust = 0) +
  fill_scale +
  labs(x = "CAZyme family", y = "Mean gene count per genome",
       tag = "A") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        plot.tag = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# ── TCDB ───────────────────────────────────────────────────────────────────────
tc_labels <- c("1" = "1: Channels/pores",
               "2" = "2: Electrochemical carriers",
               "3" = "3: Primary active",
               "4" = "4: Group translocators",
               "8" = "8: Accessory factors",
               "9" = "9: Incompletely characterized")

tc_blast <- read.table("tcdb_blast", sep = "\t", header = FALSE,
                       stringsAsFactors = FALSE, quote = "")
colnames(tc_blast) <- c("qseqid","sseqid","pident","length","mismatch",
                        "gapopen","qstart","qend","sstart","send",
                        "evalue","bitscore","qcovs")[seq_len(ncol(tc_blast))]
tc_blast <- tc_blast[tc_blast$evalue <= 1e-5 & tc_blast$pident >= 30, ]
tc_blast <- tc_blast |> group_by(qseqid) |>
  slice_max(bitscore, n = 1, with_ties = FALSE) |> ungroup()
tc_blast$tc_number <- sub(".*\\|", "", tc_blast$sseqid)
tc_blast$tc_class  <- substr(tc_blast$tc_number, 1, 1)
tc_blast <- tc_blast[tc_blast$tc_class %in% names(tc_labels), ]
tc_blast$genome <- sub("__.*$", "", tc_blast$qseqid)

tc_counts <- tc_blast |> group_by(genome, tc_class) |>
  summarise(count = n(), .groups = "drop")
tc_grid <- expand.grid(genome = unique(meta$FileName),
                       tc_class = names(tc_labels),
                       stringsAsFactors = FALSE)
tc_counts <- left_join(tc_grid, tc_counts, by = c("genome", "tc_class"))
tc_counts$count[is.na(tc_counts$count)] <- 0
tc_counts <- left_join(tc_counts, meta, by = c("genome" = "FileName"))
tc_counts <- tc_counts[!is.na(tc_counts$group), ]
tc_stats <- make_wilcox_stats(tc_counts, "tc_class")
tc_stats$tc_label <- tc_labels[tc_stats$tc_class]
tc_counts$tc_label <- tc_labels[tc_counts$tc_class]

tc_plot_df <- tc_counts |> group_by(tc_label, group) |>
  summarise(mean = mean(count), se = sd(count)/sqrt(n()), .groups = "drop")
tc_sig <- left_join(
  tc_plot_df |> filter(group == "highT") |> select(tc_label, mean),
  tc_stats |> select(tc_label, sig), by = "tc_label")

p_tc <- ggplot(tc_plot_df, aes(x = tc_label, y = mean, fill = group)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_errorbar(aes(ymin = pmax(mean - se, 0), ymax = mean + se),
                position = position_dodge(0.7), width = 0.25, linewidth = 0.5) +
  geom_text(data = tc_sig, aes(x = tc_label, y = mean + 1, label = sig),
            inherit.aes = FALSE, size = 4, vjust = 0) +
  fill_scale +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 14)) +
  labs(x = "TC class", y = "Mean transporter count per genome",
       tag = "B") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.tag = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

# ── MEROPS ─────────────────────────────────────────────────────────────────────
fam_map <- read.table("my_mer_families.txt", sep = "\t", header = FALSE,
                      stringsAsFactors = FALSE, fill = TRUE, quote = "")
colnames(fam_map) <- c("mer_id", "subfamily", "family")[seq_len(ncol(fam_map))]
fam_map <- fam_map[nchar(fam_map$mer_id) > 0 & nchar(fam_map$family) > 0,
                   c("mer_id", "family")]

mer_blast <- read.table("merops_blast", sep = "\t", header = FALSE,
                        stringsAsFactors = FALSE, quote = "")
colnames(mer_blast) <- c("qseqid","sseqid","pident","length","mismatch",
                         "gapopen","qstart","qend","sstart","send",
                         "evalue","bitscore","qcovs")[seq_len(ncol(mer_blast))]
mer_blast <- mer_blast[mer_blast$evalue <= 1e-5 & mer_blast$pident >= 30, ]
mer_blast <- mer_blast |> group_by(qseqid) |>
  slice_max(bitscore, n = 1, with_ties = FALSE) |> ungroup()
mer_blast <- left_join(mer_blast, fam_map, by = c("sseqid" = "mer_id"))
mer_blast <- mer_blast[!is.na(mer_blast$family) & nchar(mer_blast$family) > 0, ]
mer_blast$genome <- sub("__.*$", "", mer_blast$qseqid)

mer_counts <- mer_blast |> group_by(genome, family) |>
  summarise(count = n(), .groups = "drop")
mer_grid <- expand.grid(genome = unique(meta$FileName),
                        family = sort(unique(fam_map$family)),
                        stringsAsFactors = FALSE)
mer_counts <- left_join(mer_grid, mer_counts, by = c("genome", "family"))
mer_counts$count[is.na(mer_counts$count)] <- 0
mer_counts <- left_join(mer_counts, meta, by = c("genome" = "FileName"))
mer_counts <- mer_counts[!is.na(mer_counts$group), ]
mer_stats <- make_wilcox_stats(mer_counts, "family")
mer_stats <- mer_stats[!(mer_stats$mean_lowT == 0 & mer_stats$mean_highT == 0), ]
mer_counts <- mer_counts[mer_counts$family %in% mer_stats$family, ]

mer_plot_df <- mer_counts |> group_by(family, group) |>
  summarise(mean = mean(count), se = sd(count)/sqrt(n()), .groups = "drop")
mer_sig <- left_join(
  mer_plot_df |> filter(group == "highT") |> select(family, mean),
  mer_stats |> select(family, sig), by = "family")

p_mer <- ggplot(mer_plot_df, aes(x = family, y = mean, fill = group)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_errorbar(aes(ymin = pmax(mean - se, 0), ymax = mean + se),
                position = position_dodge(0.7), width = 0.25, linewidth = 0.5) +
  geom_text(data = mer_sig, aes(x = family, y = mean + 0.3, label = sig),
            inherit.aes = FALSE, size = 4, vjust = 0) +
  fill_scale +
  labs(x = "MEROPS family", y = "Mean peptidase count per genome",
       tag = "C") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
        plot.tag = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))

# ── Combine horizontally (shared legend, proportional widths) ─────────────────
combined <- (p_caz | p_tc | p_mer) +
  plot_layout(widths = c(1, 1.3, 2.4), guides = "collect") &
  theme(legend.position = "top")

ggsave("combined_family_barplots.pdf", combined, width = 18, height = 5.5)
ggsave("combined_family_barplots.png", combined, width = 18, height = 5.5, dpi = 300)
message("Done: combined_family_barplots.pdf / .png")
