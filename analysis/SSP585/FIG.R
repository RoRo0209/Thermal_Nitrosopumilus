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

## =========================================================
## FIG.R  —  GLOBAL nearshore 2×3 maps + summary stats
## LowT vs HighT, SSP585 2080-2100
## =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
})

## -----------------------------
## 0) USER SETTINGS
## -----------------------------
scenario     <- "SSP585"
period       <- "2080-2100"
mode         <- "KEEPCHL"
nearshore_km <- 200

base_low  <- thermal_path("SSP585/NP_lowT")
base_high <- thermal_path("SSP585/NP_highT")

out_dir_low  <- file.path(base_low,  "AOA_future_projection_outputs", "GLOBAL_NEARSHORE_MAPS")
out_dir_high <- file.path(base_high, "AOA_future_projection_outputs", "GLOBAL_NEARSHORE_MAPS")

out_pdf_dir <- out_dir_low   # save final PDFs here
dir.create(out_pdf_dir, showWarnings = FALSE, recursive = TRUE)

res_deg      <- 0.25
xlim         <- c(-180, 180)
ylim         <- c(-70,   80)
lon_breaks   <- seq(-180, 180, by = 60)
lat_breaks   <- seq(-60,   60, by = 30)
clip_quant   <- c(0.02, 0.98)
simplify_deg <- 0.12

## -----------------------------
## 1) HELPERS
## -----------------------------
tagify <- function(mode, scenario, period, nearshore_km) {
  paste0("nearshore", nearshore_km, "km_", mode, "_", scenario, "_", period)
}

read_three_tables <- function(out_dir, tag) {
  list(
    cur   = read_tsv(file.path(out_dir, paste0("GLOBAL_nearshore_pred_current_", tag, ".tsv")),
                     show_col_types = FALSE),
    fut   = read_tsv(file.path(out_dir, paste0("GLOBAL_nearshore_pred_future_",  tag, ".tsv")),
                     show_col_types = FALSE),
    delta = read_tsv(file.path(out_dir, paste0("GLOBAL_nearshore_pred_delta_",   tag, ".tsv")),
                     show_col_types = FALSE)
  )
}

prep_df <- function(df, zcol, xlim, ylim) {
  df %>%
    transmute(
      Longitude = as.numeric(Longitude),
      Latitude  = as.numeric(Latitude),
      z         = as.numeric(.data[[zcol]])
    ) %>%
    filter(
      is.finite(Longitude), is.finite(Latitude), is.finite(z),
      Longitude >= xlim[1], Longitude <= xlim[2],
      Latitude  >= ylim[1], Latitude  <= ylim[2]
    )
}

grid_df <- function(d, res_deg, xlim, ylim) {
  d %>%
    mutate(
      Longitude = floor((Longitude - xlim[1]) / res_deg) * res_deg + xlim[1] + res_deg / 2,
      Latitude  = floor((Latitude  - ylim[1]) / res_deg) * res_deg + ylim[1] + res_deg / 2
    ) %>%
    group_by(Longitude, Latitude) %>%
    summarise(z = mean(z, na.rm = TRUE), .groups = "drop")
}

## -----------------------------
## 2) BASE MAP (sf)
## -----------------------------
sf::sf_use_s2(FALSE)

world_sf <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf") %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  st_simplify(dTolerance = simplify_deg, preserveTopology = TRUE)

base_map <- list(
  geom_sf(data = world_sf, fill = "grey92", color = "grey55",
          linewidth = 0.2, inherit.aes = FALSE)
)

coord_ll <- coord_sf(xlim = xlim, ylim = ylim, expand = FALSE)

scale_ll <- list(
  scale_x_continuous(breaks = lon_breaks, labels = function(x) paste0(x, "°")),
  scale_y_continuous(breaks = lat_breaks, labels = function(y) paste0(y, "°"))
)

theme_pub <- theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "#D6EAF8", color = NA),
    panel.border     = element_rect(fill = NA, color = "black", linewidth = 0.4),
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 10),
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    axis.text        = element_text(size = 7)
  )

## -----------------------------
## 3) PLOT FUNCTIONS
## -----------------------------
plot_abundance <- function(df, zcol, title, lims, show_legend = TRUE) {
  d <- prep_df(df, zcol, xlim, ylim) %>% grid_df(res_deg, xlim, ylim)
  p <- ggplot(d, aes(Longitude, Latitude, fill = z)) +
    geom_tile(width = res_deg, height = res_deg) +
    base_map + coord_ll + scale_ll +
    scale_fill_viridis_c(limits = lims, oob = squish,
                         name = expression(log[10]~abundance)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_pub
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

plot_delta <- function(df, zcol, title, lims, show_legend = TRUE) {
  d <- prep_df(df, zcol, xlim, ylim) %>% grid_df(res_deg, xlim, ylim)
  p <- ggplot(d, aes(Longitude, Latitude, fill = z)) +
    geom_tile(width = res_deg, height = res_deg) +
    base_map + coord_ll + scale_ll +
    scale_fill_gradient2(limits = lims, midpoint = 0,
                         low = "#4575B4", mid = "white", high = "#D73027",
                         oob = squish,
                         name = expression(Delta~log[10])) +
    labs(title = title, x = NULL, y = NULL) +
    theme_pub
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

## -----------------------------
## 4) LOAD DATA
## -----------------------------
tag  <- tagify(mode, scenario, period, nearshore_km)
low  <- read_three_tables(out_dir_low,  tag)
high <- read_three_tables(out_dir_high, tag)

## -----------------------------
## 5) UNIFIED SCALES
## -----------------------------
z_all <- c(
  prep_df(low$cur,  "pred_current", xlim, ylim)$z,
  prep_df(low$fut,  "pred_future",  xlim, ylim)$z,
  prep_df(high$cur, "pred_current", xlim, ylim)$z,
  prep_df(high$fut, "pred_future",  xlim, ylim)$z
)
lims_abun <- quantile(z_all, probs = clip_quant, na.rm = TRUE)

z_delta <- c(
  prep_df(low$delta,  "delta_future_minus_current", xlim, ylim)$z,
  prep_df(high$delta, "delta_future_minus_current", xlim, ylim)$z
)
q_delta    <- quantile(abs(z_delta), 0.98, na.rm = TRUE)
lims_delta <- c(-q_delta, q_delta)

message("[SCALE] abundance lims: ", paste(round(lims_abun, 3), collapse=" ~ "))
message("[SCALE] delta    lims: ", paste(round(lims_delta, 3), collapse=" ~ "))

## -----------------------------
## 6) BUILD 2x3 MAP PANELS
## -----------------------------
p11 <- plot_abundance(low$cur,  "pred_current",
                      paste0("LowT  Current"),  lims_abun, show_legend = TRUE)
p12 <- plot_abundance(high$cur, "pred_current",
                      paste0("HighT Current"),  lims_abun, show_legend = FALSE)

p21 <- plot_abundance(low$fut,  "pred_future",
                      paste0("LowT  Future (", scenario, " ", period, ")"), lims_abun, show_legend = FALSE)
p22 <- plot_abundance(high$fut, "pred_future",
                      paste0("HighT Future (", scenario, " ", period, ")"), lims_abun, show_legend = FALSE)

p31 <- plot_delta(low$delta,  "delta_future_minus_current",
                  "LowT  Delta (future - current)", lims_delta, show_legend = TRUE)
p32 <- plot_delta(high$delta, "delta_future_minus_current",
                  "HighT Delta (future - current)", lims_delta, show_legend = FALSE)

combined <- (p11 | p12) / (p21 | p22) / (p31 | p32) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

out_map_pdf <- file.path(out_pdf_dir,
  paste0("GLOBAL_nearshore_2x3_", tag, ".pdf"))

ggsave(out_map_pdf, combined,
       width = 14, height = 13, dpi = 300,
       device = cairo_pdf, bg = "white")
message("Map PDF saved: ", out_map_pdf)

## -----------------------------
## 7) STATS: area-weighted total abundance
## -----------------------------
area_weight <- function(lat) cos(lat * pi / 180)

calc_total <- function(df, zcol) {
  df %>%
    transmute(lat = Latitude, v = as.numeric(.data[[zcol]])) %>%
    filter(is.finite(v)) %>%
    mutate(linear = 10^v, w = area_weight(lat)) %>%
    summarise(total = sum(linear * w, na.rm = TRUE)) %>%
    pull(total)
}

total_change <- tibble(
  group   = c("LowT", "HighT"),
  current = c(calc_total(low$cur,  "pred_current"),
              calc_total(high$cur, "pred_current")),
  future  = c(calc_total(low$fut,  "pred_future"),
              calc_total(high$fut, "pred_future"))
) %>%
  mutate(ratio = future / current,
         pct_change = (ratio - 1) * 100)

message("\n[TOTAL ABUNDANCE CHANGE]")
print(total_change)

## -----------------------------
## 8) STATS: latitudinal profile
## -----------------------------
lat_profile <- bind_rows(
  low$delta  %>% mutate(group = "LowT"),
  high$delta %>% mutate(group = "HighT")
) %>%
  transmute(group, lat = Latitude,
            delta = delta_future_minus_current,
            lat_band = floor(lat / 5) * 5) %>%
  filter(is.finite(delta)) %>%
  group_by(group, lat_band) %>%
  summarise(median_delta = median(delta, na.rm = TRUE),
            q25 = quantile(delta, 0.25, na.rm = TRUE),
            q75 = quantile(delta, 0.75, na.rm = TRUE),
            .groups = "drop")

## -----------------------------
## 9) STATS: significant change area
## -----------------------------
sig_area <- bind_rows(
  low$delta  %>% mutate(group = "LowT"),
  high$delta %>% mutate(group = "HighT")
) %>%
  transmute(group, delta = delta_future_minus_current,
            w = area_weight(Latitude),
            change_class = case_when(
              delta >  0.2 ~ "Increase",
              delta < -0.2 ~ "Decrease",
              TRUE         ~ "Stable")) %>%
  filter(is.finite(delta)) %>%
  group_by(group, change_class) %>%
  summarise(area = sum(w, na.rm = TRUE), .groups = "drop") %>%
  group_by(group) %>%
  mutate(frac = area / sum(area))

## -----------------------------
## 10) SUMMARY FIGURE A/B/C
## -----------------------------
pA <- ggplot(total_change, aes(x = group, y = ratio, fill = group)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c(LowT = "#8EC6E8", HighT = "#F4A261")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  labs(x = NULL, y = "Future / Current\n(area-weighted)",
       title = "Integrated abundance change") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

pB <- ggplot(lat_profile,
             aes(x = lat_band, y = median_delta, color = group, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.25, color = NA) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(-60, 60, 20),
                     labels = function(x) paste0(x, "°")) +
  scale_color_manual(values = c(LowT = "#8EC6E8", HighT = "#F4A261")) +
  scale_fill_manual( values = c(LowT = "#8EC6E8", HighT = "#F4A261")) +
  labs(x = "Latitude", y = expression(Delta~log[10]~abundance),
       title = "Latitudinal pattern of change",
       color = NULL, fill = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

pC <- ggplot(sig_area, aes(x = group, y = frac, fill = change_class)) +
  geom_col(width = 0.6) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c(Increase = "#D73027", Stable = "grey80",
                               Decrease = "#4575B4")) +
  labs(x = NULL, y = "Area fraction", fill = "Change",
       title = "Extent of significant change") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

summary_fig <- (pA | pB | pC) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 13)))

out_summary_pdf <- file.path(out_pdf_dir,
  paste0("SUMMARY_FIG_LowT_vs_HighT_", tag, ".pdf"))

ggsave(out_summary_pdf, summary_fig,
       width = 13, height = 4.5, dpi = 300,
       device = cairo_pdf, bg = "white")

message("Summary PDF saved: ", out_summary_pdf)
message("\nDONE ✅")
