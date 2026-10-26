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

## Supplementary figure for LowT environmental context under SSP5-8.5.
## Base R only. This mirrors the HighT SST-chlorophyll stratification logic
## used in Fig. 8f,g, but is kept as a supplementary LowT-specific analysis.

base_dir <- thermal_path("SSP585")
tag <- "nearshore200km_KEEPCHL_SSP585_2080-2100"
low_map <- file.path(
  base_dir, "NP_lowT", "AOA_future_projection_outputs", "GLOBAL_NEARSHORE_MAPS"
)
delta_file <- file.path(low_map, paste0("GLOBAL_nearshore_pred_delta_", tag, ".tsv"))

out_stem <- file.path(base_dir, "FigS_LowT_environmental_context")

sst_cut <- 22
chl_cut <- 1.5
delta_cut <- 0.2

d <- read.delim(delta_file, check.names = FALSE, stringsAsFactors = FALSE)
needed <- c("temp", "chl", "delta_future_minus_current")
if (!all(needed %in% names(d))) {
  stop("Missing required columns in ", delta_file)
}

df <- data.frame(
  future_sst = as.numeric(d$temp),
  future_chl = as.numeric(d$chl),
  delta_log10 = as.numeric(d$delta_future_minus_current)
)
df <- df[is.finite(df$future_sst) & is.finite(df$future_chl) & is.finite(df$delta_log10), ]

df$change_class <- ifelse(
  df$delta_log10 > delta_cut, "Increase",
  ifelse(df$delta_log10 < -delta_cut, "Decrease", "Stable")
)
df$change_class <- factor(df$change_class, levels = c("Increase", "Stable", "Decrease"))

sector <- function(sst, chl) {
  ifelse(
    sst >= sst_cut & chl >= chl_cut, "Warm + Productive",
    ifelse(
      sst >= sst_cut & chl < chl_cut, "Warm + Low-productivity",
      ifelse(sst < sst_cut & chl >= chl_cut, "Cool + Productive", "Cool + Low-productivity")
    )
  )
}
df$sector <- sector(df$future_sst, df$future_chl)
sector_order <- c(
  "Warm + Productive",
  "Warm + Low-productivity",
  "Cool + Productive",
  "Cool + Low-productivity"
)
df$sector <- factor(df$sector, levels = sector_order)

classes <- levels(df$change_class)
class_counts <- table(df$change_class)

quad_counts <- as.data.frame(table(df$change_class, df$sector), stringsAsFactors = FALSE)
names(quad_counts) <- c("change_class", "sector", "n")
quad_counts$total_in_class <- ave(quad_counts$n, quad_counts$change_class, FUN = sum)
quad_counts$fraction <- quad_counts$n / quad_counts$total_in_class
quad_counts$percent <- quad_counts$fraction * 100

write.table(
  quad_counts,
  file.path(base_dir, "FigS_LowT_environmental_context_sector_source.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

class_summary <- aggregate(
  cbind(delta_log10, future_sst, future_chl) ~ change_class,
  data = df,
  FUN = median
)
class_summary$n <- as.integer(class_counts[as.character(class_summary$change_class)])
class_summary$percent_of_cells <- class_summary$n / nrow(df) * 100
write.table(
  class_summary,
  file.path(base_dir, "FigS_LowT_environmental_context_class_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

sector_cols <- c(
  "Warm + Productive" = "#D73027",
  "Warm + Low-productivity" = "#FDAE61",
  "Cool + Productive" = "#74ADD1",
  "Cool + Low-productivity" = "#ABD9E9"
)
class_cols <- c("Increase" = "#C95F42", "Stable" = "#808080", "Decrease" = "#3E76A8")

alpha_col <- function(col, alpha = 0.18) {
  rgb_col <- col2rgb(col) / 255
  rgb(rgb_col[1], rgb_col[2], rgb_col[3], alpha = alpha)
}

draw_sector_labels <- function(xlim, ylim) {
  text(xlim[2] - 1.2, ylim[2] * 0.92, "Warm +\nProductive", adj = c(1, 1),
       cex = 0.47, col = "#333333", font = 3)
  text(xlim[2] - 1.2, ylim[1] + 0.18, "Warm +\nLow-productivity", adj = c(1, 0),
       cex = 0.47, col = "#333333", font = 3)
  text(xlim[1] + 1.2, ylim[2] * 0.92, "Cool +\nProductive", adj = c(0, 1),
       cex = 0.47, col = "#333333", font = 3)
  text(xlim[1] + 1.2, ylim[1] + 0.18, "Cool +\nLow-productivity", adj = c(0, 0),
       cex = 0.47, col = "#333333", font = 3)
}

draw_figure <- function() {
  old <- par(
    oma = c(0.2, 0.2, 0.2, 0.2),
    family = "sans",
    lend = "round",
    las = 1
  )
  on.exit(par(old))
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE), heights = c(2.15, 1.05))

  xlim <- c(-2, 38)
  ylim <- c(0, as.numeric(quantile(df$future_chl, 0.98, na.rm = TRUE)))

  for (cl in classes) {
    sub <- df[df$change_class == cl, ]
    par(mar = c(3.0, 3.15, 1.55, 0.45), mgp = c(1.75, 0.46, 0), tcl = -0.18,
        cex.axis = 0.58, cex.lab = 0.66)
    plot(
      sub$future_sst,
      pmin(sub$future_chl, ylim[2]),
      xlim = xlim,
      ylim = ylim,
      pch = 16,
      cex = 0.18,
      col = alpha_col(class_cols[[as.character(cl)]], 0.18),
      axes = FALSE,
      xlab = "",
      ylab = "",
      bty = "n"
    )
    abline(v = sst_cut, lty = 2, lwd = 0.8, col = "#333333")
    abline(h = chl_cut, lty = 2, lwd = 0.8, col = "#333333")
    draw_sector_labels(xlim, ylim)
    axis(1, at = seq(0, 40, 10))
    axis(2, at = seq(0, floor(ylim[2]), 1))
    box(bty = "o", lwd = 0.85, col = "black")
    mtext("Future SST (°C)", side = 1, line = 1.58, cex = 0.62)
    if (cl == "Increase") {
      mtext(expression("Future Chl (mg m"^{-3}*")"), side = 2, line = 2.05, las = 0, cex = 0.62)
    }
    title(
      main = paste0(as.character(cl), " (n = ", format(class_counts[as.character(cl)], big.mark = ","), ")"),
      cex.main = 0.68,
      font.main = 2,
      line = 0.45
    )
  }

  par(mar = c(2.95, 3.15, 0.35, 2.7), mgp = c(1.75, 0.46, 0), tcl = -0.18,
      cex.axis = 0.58, cex.lab = 0.66)
  plot(
    NA,
    xlim = c(0.4, 3.6),
    ylim = c(0, 100),
    axes = FALSE,
    xlab = "",
    ylab = "",
    bty = "n"
  )
  bar_x <- 1:3
  bar_w <- 0.48
  for (i in seq_along(classes)) {
    cl <- classes[i]
    sub <- quad_counts[quad_counts$change_class == cl, ]
    sub <- sub[match(sector_order, sub$sector), ]
    y0 <- 0
    for (j in seq_len(nrow(sub))) {
      y1 <- y0 + sub$percent[j]
      rect(bar_x[i] - bar_w / 2, y0, bar_x[i] + bar_w / 2, y1,
           col = sector_cols[as.character(sub$sector[j])], border = "white", lwd = 0.25)
      if (sub$percent[j] >= 7) {
        text(bar_x[i], (y0 + y1) / 2, paste0(round(sub$percent[j]), "%"),
             cex = 0.54, col = ifelse(j == 1, "white", "#222222"))
      }
      y0 <- y1
    }
  }
  axis(1, at = bar_x, labels = classes)
  axis(2, at = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"))
  box(bty = "o", lwd = 0.85, col = "black")
  mtext("LowT abundance change", side = 1, line = 1.62, cex = 0.62)
  mtext("Fraction of grid cells", side = 2, line = 2.05, las = 0, cex = 0.62)
  legend(
    "right",
    inset = c(-0.22, 0),
    legend = names(sector_cols),
    fill = sector_cols,
    border = NA,
    bty = "n",
    cex = 0.55,
    title = "Environmental sector",
    xpd = TRUE
  )
}

pdf(paste0(out_stem, ".pdf"), width = 7.2, height = 5.2,
    family = "Helvetica", useDingbats = FALSE)
draw_figure()
dev.off()

svg(paste0(out_stem, ".svg"), width = 7.2, height = 5.2,
    family = "sans", pointsize = 10)
draw_figure()
dev.off()

tiff(paste0(out_stem, ".tiff"), width = 7.2, height = 5.2,
     units = "in", res = 600, compression = "lzw", family = "sans")
draw_figure()
dev.off()

png(paste0(out_stem, ".png"), width = 7.2, height = 5.2,
    units = "in", res = 300, family = "sans")
draw_figure()
dev.off()

cat("Wrote:", out_stem, "{pdf,svg,tiff,png}\n")
print(quad_counts)
