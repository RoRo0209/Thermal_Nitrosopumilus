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

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## -----------------------------
## 0) Paths
## -----------------------------
indir  <- thermal_path("genome_character/two_group")
infile <- file.path(indir, "two_group.xlsx")

out_png <- file.path(indir, "two_group_equalwidth_fixedorder_with_mu.png")
out_pdf <- file.path(indir, "two_group_equalwidth_fixedorder_with_mu.pdf")

stopifnot(file.exists(infile))

## -----------------------------
## 1) Read data
## -----------------------------
df <- read_excel(infile) %>% as.data.frame()

## -----------------------------
## 2) group filter + order
## -----------------------------
df <- df %>%
  mutate(group = trimws(as.character(group))) %>%
  filter(group %in% c("lowT", "highT")) %>%
  mutate(group = factor(group, levels = c("lowT", "highT")))

cat("\n[INFO] group counts:\n")
print(table(df$group))

## -----------------------------
## 3) numeric columns (genome traits)
## -----------------------------
df <- df %>%
  mutate(
    OGT            = as.numeric(predOGT),
    GC_pct         = as.numeric(`GC-content`),
    Genome_Mbp     = as.numeric(Genome_Size) / 1e6,   # 如果本来就是 Mbp，请删掉 /1e6
    CAZymes        = as.numeric(dbcan),
    Peptidases     = as.numeric(merops),
    Coding_Density = as.numeric(Coding_Density),
    TCDB           = as.numeric(tcdb)
  )

## coding density 如果是 0–1，转 %
coding_unit <- "%"
if ("Coding_Density" %in% names(df) && max(df$Coding_Density, na.rm = TRUE) <= 1.2) {
  df <- df %>% mutate(Coding_Density = Coding_Density * 100)
}

## -----------------------------
## 4) gRodon columns (SAFE & robust)
## -----------------------------

# 安全取列：列不存在就返回 NA_real_ 向量（长度 = nrow(df)）
get_num <- function(dat, col) {
  if (col %in% names(dat)) return(as.numeric(dat[[col]]))
  rep(NA_real_, nrow(dat))
}

d1  <- get_num(df, "d_hours")
d2  <- get_num(df, "doubling_time")
d3  <- get_num(df, "d")
mu0 <- get_num(df, "mu_per_hour")

# 先合并出 doubling time（小时）
d_hours_vec <- d1
d_hours_vec[is.na(d_hours_vec)] <- d2[is.na(d_hours_vec)]
d_hours_vec[is.na(d_hours_vec)] <- d3[is.na(d_hours_vec)]

# 再合并 mu：优先用已有的 mu_per_hour；否则用 ln2/d
mu_vec <- mu0
mu_vec[is.na(mu_vec) & !is.na(d_hours_vec) & d_hours_vec > 0] <- log(2) / d_hours_vec[is.na(mu_vec) & !is.na(d_hours_vec) & d_hours_vec > 0]
if (!any(is.finite(mu_vec))) {
  stop("Growth-rate panel requires finite mu_per_hour or positive doubling-time input.")
}

df <- df %>%
  mutate(
    d_hours = d_hours_vec,
    mu_per_hour = mu_vec
  )

cat("\n[INFO] d_hours NA fraction:", mean(is.na(df$d_hours)), "\n")
cat("[INFO] mu_per_hour NA fraction:", mean(is.na(df$mu_per_hour)), "\n")


## sanity: 确认 mu_per_hour 不是全 NA
cat("\n[INFO] mu_per_hour NA fraction:\n")
print(mean(is.na(df$mu_per_hour)))

## -----------------------------
## 5) colors
## -----------------------------
col_fill <- c("lowT"= "#8EC6E8", "highT"="#F4A261")
col_pt   <- c("lowT"= "#8EC6E8", "highT"="#F4A261")

p_to_star <- function(p){
  if (is.na(p)) return("ns")
  if (p < 0.0005) return("***")
  if (p < 0.005)  return("**")
  if (p < 0.05)   return("*")
  return("ns")
}

## -----------------------------
## 6) panel maker
## -----------------------------
make_panel <- function(dat, y, ylab, title, show_legend = FALSE){
  
  stopifnot(nlevels(dat$group) == 2)
  
  yvec <- dat[[y]]
  wt <- wilcox.test(yvec ~ dat$group, exact = FALSE)
  star <- p_to_star(wt$p.value)
  
  ymax <- max(yvec, na.rm = TRUE)
  ymin <- min(yvec, na.rm = TRUE)
  yrng <- ymax - ymin
  if (!is.finite(yrng) || yrng == 0) yrng <- 1
  y_star <- ymax + 0.12 * yrng
  
  p <- ggplot(dat, aes(x = group, y = .data[[y]], fill = group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA,
                 color = "black", linewidth = 0.7) +
    geom_jitter(aes(color = group),
                width = 0.12, size = 2.2, alpha = 0.95) +
    annotate("segment", x = 1, xend = 2,
             y = y_star, yend = y_star,
             linewidth = 1.1) +
    annotate("text", x = 1.5,
             y = y_star + 0.03 * yrng,
             label = star, size = 6, fontface = "bold") +
    scale_fill_manual(values = col_fill) +
    scale_color_manual(values = col_pt) +
    labs(title = title, x = NULL, y = ylab, fill = "group") +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      plot.title  = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 6)),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1.0),
      plot.margin = margin(6, 6, 6, 6)
    )
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 13)
    )
  }
  p
}

## -----------------------------
## 7) build panels (顺序严格固定 2x4)
## -----------------------------
p1 <- make_panel(df, "OGT", "Degrees Celsius", "Optimal growth temperature")
p2 <- make_panel(df, "GC_pct", "GC %", "GC-content")
p3 <- make_panel(df, "CAZymes", "Genes per genome", "CAZymes")
p4 <- make_panel(df, "Peptidases", "Genes per genome", "Peptidases", show_legend = TRUE)

p5 <- make_panel(df, "Genome_Mbp", "Megabase pairs", "Genome size")
p6 <- make_panel(df, "Coding_Density", paste0("Coding density (", coding_unit, ")"), "Coding density")
p7 <- make_panel(df, "TCDB", "Genes per genome", "TCDB")

## NEW: p8 = mu_per_hour
## 如果 mu_per_hour 全 NA，这里会在 wilcox.test 报错；提前 stop 更清晰
if (all(is.na(df$mu_per_hour))) {
  stop("mu_per_hour is all NA. Please check your gRodon columns (d/d_hours/mu_per_hour).")
}
p8 <- make_panel(df, "mu_per_hour", expression(paste(mu, " (h"^{-1}, ")")), "Max growth rate (gRodon)")

final <- (p1 + p2 + p3 + p4) /
  (p5 + p6 + p7 + p8) +
  plot_layout(ncol = 4)

final

## -----------------------------
## 8) save
## -----------------------------
ggsave(out_png, final, width = 18.5, height = 8.2, dpi = 300)
ggsave(out_pdf, final, width = 24, height = 8.2)

message("Done!\n", out_png, "\n", out_pdf)
