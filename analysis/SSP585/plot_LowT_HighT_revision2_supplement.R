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
setwd(thermal_path("SSP585"))

## Supplementary revision figure for paired LowT and HighT SSP5-8.5 responses.
## Panel a follows the reviewer suggestion: x = latitude-sorted grid cells,
## y = projected abundance change. Panel b quantifies the direction pairs
## across the original grid cells.

base_dir <- thermal_path("SSP585")
in_file <- file.path(base_dir, "LowT_HighT_latitude_sorted_change.tsv")
out_stem <- file.path(base_dir, "Fig_revision2_LowT_HighT_spatial_divergence_supplement")

dat <- read.delim(in_file, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c(
  "latitude_sorted_index", "Latitude", "LowT_delta_log10",
  "HighT_delta_log10", "direction_pair"
)
stopifnot(all(needed %in% names(dat)))

dat <- dat[order(dat$latitude_sorted_index), ]
dat$sorted_position <- seq_len(nrow(dat)) / nrow(dat) * 100

pair_levels <- c(
  "LowT_down_HighT_up",
  "LowT_up_HighT_down",
  "both_increase",
  "both_decrease"
)
pair_labels <- c(
  "LowT decrease / HighT increase",
  "LowT increase / HighT decrease",
  "Both increase",
  "Both decrease"
)
pair_cols <- c("#D97745", "#4F83B5", "#77A88D", "#8B8B8B")
pair_counts <- table(factor(dat$direction_pair, levels = pair_levels))
pair_pct <- as.numeric(pair_counts) / nrow(dat) * 100
opposite_pct <- sum(pair_pct[1:2])
same_pct <- sum(pair_pct[3:4])

## A lightly binned median is overlaid only as a visual guide; counts in panel b
## are calculated from the original 42,258 grid cells.
n_bins <- 240
dat$plot_bin <- ceiling(seq_len(nrow(dat)) / nrow(dat) * n_bins)
trend <- aggregate(
  cbind(Latitude, sorted_position, LowT_delta_log10, HighT_delta_log10) ~ plot_bin,
  data = dat,
  FUN = median
)

write.table(
  data.frame(
    response_pair = pair_labels,
    grid_cells = as.numeric(pair_counts),
    percent = pair_pct
  ),
  file.path(base_dir, "Fig_revision2_LowT_HighT_spatial_divergence_supplement_pair_source.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

col_low <- "#3E93AA"
col_high <- "#D8763F"
col_zero <- "#707070"
col_lat <- "#A0A0A0"

alpha_col <- function(col, alpha = 0.16) {
  rgb_col <- col2rgb(col) / 255
  rgb(rgb_col[1], rgb_col[2], rgb_col[3], alpha = alpha)
}

draw_figure <- function() {
  layout(matrix(c(1, 2), ncol = 1), heights = c(3.6, 1.35))
  old <- par(
    mar = c(3.6, 4.6, 1.0, 1.0),
    oma = c(0.2, 0.2, 0.1, 0.1),
    mgp = c(2.35, 0.62, 0),
    tcl = -0.22,
    family = "sans",
    las = 1,
    cex.axis = 0.78,
    cex.lab = 0.92,
    lend = "round"
  )
  on.exit(par(old))

  ## Panel a: paired latitude-sorted grid-cell changes.
  ylim <- range(c(dat$LowT_delta_log10, dat$HighT_delta_log10), finite = TRUE)
  ylim <- c(max(-0.45, floor(ylim[1] * 10) / 10), min(1.20, ceiling(ylim[2] * 10) / 10))
  plot(
    dat$sorted_position, dat$LowT_delta_log10,
    type = "n",
    xlim = c(0, 100),
    ylim = ylim,
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    xlab = "",
    ylab = expression(Delta * log[10] ~ abundance)
  )
  abline(h = 0, col = col_zero, lty = 2, lwd = 0.9)
  points(dat$sorted_position, dat$LowT_delta_log10,
         pch = 16, cex = 0.12, col = alpha_col(col_low, 0.10))
  points(dat$sorted_position, dat$HighT_delta_log10,
         pch = 16, cex = 0.12, col = alpha_col(col_high, 0.10))
  lines(trend$sorted_position, trend$LowT_delta_log10, col = col_low, lwd = 1.8)
  lines(trend$sorted_position, trend$HighT_delta_log10, col = col_high, lwd = 1.8)
  axis(2)
  axis(1, at = seq(0, 100, 20), labels = seq(0, 100, 20))
  mtext("Latitude-sorted coastal grid cells (%)", side = 1, line = 2.0, cex = 0.88)
  box(bty = "l", lwd = 0.8)
  legend(
    "topleft",
    legend = c("LowT", "HighT"),
    col = c(col_low, col_high),
    lwd = 2.0,
    bty = "n",
    horiz = TRUE,
    cex = 0.82,
    seg.len = 2.3,
    inset = c(0.005, -0.005)
  )
  text(1.2, ylim[2] - diff(ylim) * 0.07, "a", font = 2, cex = 1.05, adj = c(0, 0.5))
  text(
    98.5, ylim[1] + diff(ylim) * 0.08,
    "Lines show binned medians; direction classes use original grid cells",
    cex = 0.55, adj = c(1, 0.5), col = "#555555"
  )

  lat_breaks <- c(-60, -30, 0, 30, 60)
  lat_pos <- sapply(lat_breaks, function(v) {
    trend$sorted_position[which.min(abs(trend$Latitude - v))]
  })
  axis(
    3,
    at = lat_pos,
    labels = paste0(lat_breaks, "\u00b0"),
    tick = FALSE,
    line = -0.25,
    cex.axis = 0.70,
    col.axis = col_lat
  )

  ## Panel b: original grid-cell direction-pair composition.
  par(mar = c(3.0, 4.6, 0.1, 1.0))
  bar_x <- cumsum(c(0, pair_pct))
  short_labels <- c(
    "LowT down /\nHighT up",
    "LowT up /\nHighT down",
    "Both\nincrease",
    "Both\ndecrease"
  )
  plot(
    NA,
    xlim = c(0, 100),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    xlab = "",
    ylab = "",
    bty = "n"
  )
  for (i in seq_along(pair_pct)) {
    rect(bar_x[i], 0.48, bar_x[i + 1], 0.80,
         col = pair_cols[i], border = "white", lwd = 1.1)
    text(
      (bar_x[i] + bar_x[i + 1]) / 2,
      0.64,
      sprintf("%.1f%%", pair_pct[i]),
      cex = 0.77,
      col = if (i == 4) "white" else "#242424"
    )
  }
  axis(1, at = seq(0, 100, 20), labels = seq(0, 100, 20))
  mtext("Grid cells (%)", side = 1, line = 2.0, cex = 0.88)
  text(1.2, 0.94, "b", font = 2, cex = 1.05, adj = c(0, 0.5))
  text(
    98.5, 0.94,
    sprintf("Opposite directions %.1f%%; same directions %.1f%%", opposite_pct, same_pct),
    cex = 0.62, adj = c(1, 0.5), col = "#444444"
  )
  text(
    (bar_x[-length(bar_x)] + bar_x[-1]) / 2,
    0.27,
    short_labels,
    cex = 0.61,
    col = "#444444"
  )
}

pdf(paste0(out_stem, ".pdf"), width = 7.2, height = 5.6,
    family = "Helvetica", useDingbats = FALSE)
draw_figure()
dev.off()

svg(paste0(out_stem, ".svg"), width = 7.2, height = 5.6,
    family = "sans", pointsize = 10)
draw_figure()
dev.off()

tiff(paste0(out_stem, ".tiff"), width = 7.2, height = 5.6,
     units = "in", res = 600, compression = "lzw", family = "sans")
draw_figure()
dev.off()

png(paste0(out_stem, ".png"), width = 7.2, height = 5.6,
    units = "in", res = 300, family = "sans")
draw_figure()
dev.off()

cat("Wrote:", out_stem, "{pdf,svg,tiff,png}\n")
cat(sprintf("Opposite: %.2f%% | Same: %.2f%%\n", opposite_pct, same_pct))
