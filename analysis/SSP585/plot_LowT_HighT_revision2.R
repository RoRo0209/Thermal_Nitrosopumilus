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

## Revision 2 figure: paired LowT and HighT SSP5-8.5 responses.
## Base R only. The top panel summarizes responses in latitude bands. The
## bottom panel quantifies all 42,258 original grid-cell direction pairs.

base_dir <- thermal_path("SSP585")
in_file <- file.path(base_dir, "LowT_HighT_latitude_sorted_change.tsv")
out_stem <- file.path(base_dir, "Fig_revision2_LowT_HighT_spatial_divergence_revised")

dat <- read.delim(in_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(
  all(c("latitude_sorted_index", "Latitude", "LowT_delta_log10",
        "HighT_delta_log10", "direction_pair") %in% names(dat))
)

## Summarize grid cells in 5-degree latitude bands so the x axis is directly
## interpretable as geographic latitude.
lat_bin_width <- 5
dat$plot_bin <- floor(dat$Latitude / lat_bin_width) * lat_bin_width + lat_bin_width / 2
trend <- aggregate(
  cbind(LowT_delta_log10, HighT_delta_log10) ~ plot_bin,
  data = dat,
  FUN = median
)
names(trend)[1] <- "Latitude"
trend <- trend[order(trend$Latitude), ]

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

source_trend <- trend[, c("Latitude", "LowT_delta_log10", "HighT_delta_log10")]
write.table(
  source_trend,
  file.path(base_dir, "Fig_revision2_LowT_HighT_spatial_divergence_trend_source.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(
    response_pair = pair_labels,
    grid_cells = as.numeric(pair_counts),
    percent = pair_pct
  ),
  file.path(base_dir, "Fig_revision2_LowT_HighT_spatial_divergence_pair_source.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

col_low <- "#3E93AA"
col_high <- "#D8763F"
col_zero <- "#707070"
col_lat <- "#A0A0A0"

draw_figure <- function() {
  layout(matrix(c(1, 2), ncol = 1), heights = c(3.5, 1.35))
  old <- par(
    mar = c(3.5, 4.4, 1.0, 1.0),
    oma = c(0.2, 0.2, 0.1, 0.1),
    mgp = c(2.25, 0.62, 0),
    tcl = -0.22,
    family = "sans",
    las = 1,
    cex.axis = 0.80,
    cex.lab = 0.92,
    lend = "round"
  )
  on.exit(par(old))

  ## Panel a: latitude-band spatial trend
  x <- trend$Latitude
  ylim <- range(c(trend$LowT_delta_log10, trend$HighT_delta_log10, 0), na.rm = TRUE)
  ylim <- ylim + c(-0.10, 0.10) * diff(ylim)
  plot(
    x, trend$LowT_delta_log10,
    type = "n",
    xlim = c(-80, 85),
    ylim = ylim,
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    xlab = "",
    ylab = expression(Delta * log[10] ~ abundance)
  )
  abline(h = 0, col = col_zero, lty = 2, lwd = 0.9)
  polygon(
    c(x, rev(x)),
    c(trend$LowT_delta_log10, rev(trend$HighT_delta_log10)),
    col = "#E9EEF1",
    border = NA
  )
  lines(x, trend$LowT_delta_log10, col = col_low, lwd = 1.8)
  lines(x, trend$HighT_delta_log10, col = col_high, lwd = 1.8)
  points(x, trend$LowT_delta_log10, col = col_low, pch = 16, cex = 0.48)
  points(x, trend$HighT_delta_log10, col = col_high, pch = 16, cex = 0.48)
  axis(2)
  axis(1, at = seq(-60, 60, 30),
       labels = c("60°S", "30°S", "0°", "30°N", "60°N"))
  mtext("Latitude", side = 1, line = 2.0,
        cex = 0.88)
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
  text(-77, ylim[2] - 0.03 * diff(ylim), "a", font = 2, cex = 1.05,
       adj = c(0, 0.5))

  ## Panel b: direction-pair composition
  par(mar = c(3.0, 4.4, 1.25, 1.0))
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
    ylim = c(0, 1.18),
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
  segments(0.8, 0.96, opposite_pct - 0.8, 0.96, lwd = 1.0, col = "#333333")
  segments(c(0.8, opposite_pct - 0.8), 0.92,
           c(0.8, opposite_pct - 0.8), 0.96, lwd = 1.0, col = "#333333")
  text(opposite_pct / 2, 1.06,
       sprintf("Opposite directions  %.1f%%", opposite_pct),
       cex = 0.76, font = 2, col = "#333333")
  segments(opposite_pct + 0.8, 0.96, 99.2, 0.96, lwd = 1.0, col = "#333333")
  segments(c(opposite_pct + 0.8, 99.2), 0.92,
           c(opposite_pct + 0.8, 99.2), 0.96, lwd = 1.0, col = "#333333")
  text(opposite_pct + same_pct / 2, 1.06,
       sprintf("Same direction  %.1f%%", same_pct),
       cex = 0.76, font = 2, col = "#333333")
  axis(1, at = seq(0, 100, 20), labels = seq(0, 100, 20))
  mtext("Grid cells (%)", side = 1, line = 2.0, cex = 0.88)
  text(1.2, 1.10, "b", font = 2, cex = 1.05, adj = c(0, 0.5))
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
