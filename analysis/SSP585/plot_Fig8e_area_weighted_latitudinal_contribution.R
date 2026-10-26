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

## Replacement panel for Fig. 8e.
## Base R only. The panel shows area-weighted latitudinal contribution to
## global abundance change, rather than median Delta log10 abundance.

scenario <- "SSP585"
period <- "2080-2100"
mode <- "KEEPCHL"
nearshore_km <- 200

base_dir <- thermal_path("SSP585")
base_low <- file.path(base_dir, "NP_lowT")
base_high <- file.path(base_dir, "NP_highT")
tag <- paste0("nearshore", nearshore_km, "km_", mode, "_", scenario, "_", period)

read_delta <- function(base, group) {
  file <- file.path(
    base,
    "AOA_future_projection_outputs",
    "GLOBAL_NEARSHORE_MAPS",
    paste0("GLOBAL_nearshore_pred_delta_", tag, ".tsv")
  )
  d <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  needed <- c("Longitude", "Latitude", "pred_current", "pred_future")
  if (!all(needed %in% names(d))) {
    stop("Missing columns in ", file)
  }
  out <- data.frame(
    group = group,
    Longitude = as.numeric(d$Longitude),
    Latitude = as.numeric(d$Latitude),
    pred_current = as.numeric(d$pred_current),
    pred_future = as.numeric(d$pred_future)
  )
  out[is.finite(out$Longitude) & is.finite(out$Latitude) &
        is.finite(out$pred_current) & is.finite(out$pred_future), ]
}

area_weight <- function(lat) cos(lat * pi / 180)

dat <- rbind(
  read_delta(base_low, "LowT"),
  read_delta(base_high, "HighT")
)
dat$lat_band <- floor(dat$Latitude / 5) * 5
dat$current_linear <- 10^dat$pred_current
dat$future_linear <- 10^dat$pred_future
dat$weighted_delta <- (dat$future_linear - dat$current_linear) * area_weight(dat$Latitude)
dat$current_weighted <- dat$current_linear * area_weight(dat$Latitude)
dat$future_weighted <- dat$future_linear * area_weight(dat$Latitude)

groups <- c("LowT", "HighT")
lat_bands <- sort(unique(dat$lat_band))

rows <- list()
for (g in groups) {
  dg <- dat[dat$group == g, ]
  total_current <- sum(dg$current_weighted, na.rm = TRUE)
  total_future <- sum(dg$future_weighted, na.rm = TRUE)
  for (lb in lat_bands) {
    x <- dg[dg$lat_band == lb, ]
    if (nrow(x) == 0) next
    weighted_delta_sum <- sum(x$weighted_delta, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      group = g,
      lat_band = lb,
      lat_mid = lb + 2.5,
      weighted_delta_sum = weighted_delta_sum,
      current_weighted_sum = sum(x$current_weighted, na.rm = TRUE),
      future_weighted_sum = sum(x$future_weighted, na.rm = TRUE),
      n_grid = nrow(x),
      total_current_weighted = total_current,
      contribution_pct = weighted_delta_sum / total_current * 100,
      future_current_ratio = total_future / total_current
    )
  }
}
lat_contrib <- do.call(rbind, rows)

write.table(
  lat_contrib,
  file.path(base_dir, "Fig8e_area_weighted_latitudinal_contribution_source.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary_rows <- list()
for (g in groups) {
  x <- lat_contrib[lat_contrib$group == g, ]
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    group = g,
    net_weighted_delta = sum(x$weighted_delta_sum, na.rm = TRUE),
    positive_contribution_pct = sum(x$contribution_pct[x$contribution_pct > 0], na.rm = TRUE),
    negative_contribution_pct = sum(x$contribution_pct[x$contribution_pct < 0], na.rm = TRUE),
    net_contribution_pct = sum(x$contribution_pct, na.rm = TRUE),
    future_current_ratio = unique(round(x$future_current_ratio, 4))
  )
}
summary_tbl <- do.call(rbind, summary_rows)
write.table(
  summary_tbl,
  file.path(base_dir, "Fig8e_area_weighted_latitudinal_contribution_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

col_low <- "#3E93AA"
col_high <- "#D8763F"
col_zero <- "#707070"

draw_panel <- function() {
  old <- par(
    mar = c(2.75, 3.25, 0.35, 0.55),
    mgp = c(1.70, 0.42, 0),
    tcl = -0.18,
    family = "sans",
    las = 1,
    cex.axis = 0.60,
    cex.lab = 0.66,
    lend = "round"
  )
  on.exit(par(old))

  ymax <- max(abs(lat_contrib$contribution_pct), na.rm = TRUE)
  ylim <- c(-ceiling(ymax / 5) * 5, ceiling(ymax / 5) * 5)

  plot(
    NA,
    xlim = c(-70, 80),
    ylim = ylim,
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    xlab = "",
    ylab = "",
    bty = "n"
  )
  abline(h = seq(ylim[1], ylim[2], by = 5), lwd = 0.28, col = "#E6E6E6")
  abline(h = 0, lty = 2, lwd = 0.75, col = col_zero)

  low <- lat_contrib[lat_contrib$group == "LowT", ]
  high <- lat_contrib[lat_contrib$group == "HighT", ]
  bar_w <- 1.65
  rect(low$lat_mid - bar_w - 0.25, 0, low$lat_mid - 0.25, low$contribution_pct,
       col = adjustcolor(col_low, alpha.f = 0.82), border = "white", lwd = 0.15)
  rect(high$lat_mid + 0.25, 0, high$lat_mid + bar_w + 0.25, high$contribution_pct,
       col = adjustcolor(col_high, alpha.f = 0.82), border = "white", lwd = 0.15)

  axis(1, at = seq(-60, 60, 20), labels = paste0(seq(-60, 60, 20), "\u00b0"))
  axis(2, at = seq(ylim[1], ylim[2], by = 5))
  mtext("Latitude", side = 1, line = 1.55, cex = 0.66)
  mtext("Contribution to global\nabundance change (%)", side = 2, line = 2.05, las = 0, cex = 0.66)
  box(bty = "o", lwd = 0.85, col = "black")
  legend(
    "topleft",
    legend = c("LowT", "HighT"),
    fill = c(adjustcolor(col_low, alpha.f = 0.78), adjustcolor(col_high, alpha.f = 0.78)),
    border = NA,
    bty = "n",
    cex = 0.58,
    inset = c(0.015, 0.015),
    y.intersp = 0.88
  )
}

out_stem <- file.path(base_dir, "Fig8e_area_weighted_latitudinal_contribution")

panel_width <- 3.25
panel_height <- 3.0

pdf(paste0(out_stem, ".pdf"), width = panel_width, height = panel_height,
    family = "Helvetica", useDingbats = FALSE)
draw_panel()
dev.off()

svg(paste0(out_stem, ".svg"), width = panel_width, height = panel_height,
    family = "sans", pointsize = 10)
draw_panel()
dev.off()

tiff(paste0(out_stem, ".tiff"), width = panel_width, height = panel_height,
     units = "in", res = 600, compression = "lzw", family = "sans")
draw_panel()
dev.off()

png(paste0(out_stem, ".png"), width = panel_width, height = panel_height,
    units = "in", res = 300, family = "sans")
draw_panel()
dev.off()

cat("Wrote:", out_stem, "{pdf,svg,tiff,png}\n")
print(summary_tbl)
