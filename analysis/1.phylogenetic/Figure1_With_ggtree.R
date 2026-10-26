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
## Figure 1: 使用ggtree的增强版本
##
## 需要安装ggtree:
## if (!require("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install("ggtree")
## =========================================================

# 检查ggtree是否可用
if (!requireNamespace("ggtree", quietly = TRUE)) {
  cat("错误: ggtree包未安装\n")
  cat("请在Windows R中运行以下命令安装:\n")
  cat("  if (!require('BiocManager', quietly = TRUE))\n")
  cat("      install.packages('BiocManager')\n")
  cat("  BiocManager::install('ggtree')\n")
  cat("\n或使用 Figure1_Complete_Solution.R (不需要ggtree)\n")
  quit(status = 1)
}

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(patchwork)
})

setwd(thermal_path("1.phylogenetic"))

## =========================================================
## 1. 读取树文件
## =========================================================

cat("读取树文件...\n")

if (!file.exists("MCMCtree_cleaned.tre")) {
  nexus_lines <- readLines("MCMCtree_FigTree.tre")
  tree_line <- grep("UTREE", nexus_lines, value = TRUE)
  newick_str <- sub(".*= ", "", tree_line)
  newick_str <- sub(";$", "", newick_str)
  newick_clean <- gsub("\\[&[^]]+\\]", "", newick_str)
  writeLines(paste0(newick_clean, ";"), "MCMCtree_cleaned.tre")
}

tree_time <- read.tree("MCMCtree_cleaned.tre")
cat(sprintf("✓ 树文件: %d tips\n", length(tree_time$tip.label)))

## =========================================================
## 2. 读取元数据
## =========================================================

cat("读取元数据...\n")
meta <- read.xlsx("chiplot-tree.xlsx", sheet = 1)
names(meta)[names(meta) == "predOGT.(C)"] <- "OGT"

outgroup_tips <- setdiff(tree_time$tip.label, meta$FileName)

if (!is.rooted(tree_time)) {
  tree_time <- root(tree_time, outgroup = outgroup_tips, resolve.root = TRUE)
}

all_tips <- data.frame(
  FileName = tree_time$tip.label,
  stringsAsFactors = FALSE
) %>%
  left_join(meta %>% select(FileName, OGT, group), by = "FileName") %>%
  mutate(
    Clade = case_when(
      FileName %in% outgroup_tips ~ "Outgroup",
      group == "lowT" ~ "Cold-adapted",
      group == "highT" ~ "Warm-adapted",
      TRUE ~ "Other"
    )
  )

## =========================================================
## 3. 读取地质年代数据
## =========================================================

cat("读取地质年代...\n")
time2_data <- read.xlsx("TIME2.xlsx", sheet = 1)

## =========================================================
## 4. 配色
## =========================================================

isme_colors <- list(
  cold = "#3B9AB2",
  warm = "#E8601C",
  other = "#CCCCCC",
  outgroup = "#4A4A4A",
  accent = "#F21A00",
  neutral = "#2C3E50"
)

tip_colors <- c(
  "Cold-adapted" = isme_colors$cold,
  "Warm-adapted" = isme_colors$warm,
  "Other" = isme_colors$other,
  "Outgroup" = isme_colors$outgroup
)

## =========================================================
## 5. 使用ggtree绘制
## =========================================================

cat("生成ggtree图...\n")

# 基础树
p <- ggtree(tree_time, layout = "rectangular", size = 0.35) +
  theme_tree2(base_size = 8)

# 添加tip颜色
p <- p %<+% all_tips +
  geom_tippoint(aes(color = Clade), size = 1.5, alpha = 0.85) +
  scale_color_manual(values = tip_colors, name = "Ecotype", na.value = isme_colors$other)

# 反转x轴 (时间从右到左)
p <- revts(p)

# 添加地质年代背景
max_time <- max(node.depth.edgelength(tree_time))

geo_colors <- c("#FFF4E6", "#E6F2FF", "#F0F0F0", "#FFE8CC")

for (i in 1:nrow(time2_data)) {
  period <- time2_data[i, ]
  if (period$start * 100 <= max_time) {
    p <- p +
      annotate("rect",
               xmin = period$end * 100,
               xmax = min(period$start * 100, max_time),
               ymin = -Inf, ymax = Inf,
               fill = geo_colors[i], alpha = 0.15)
  }
}

# 添加年代标签
for (i in 1:nrow(time2_data)) {
  period <- time2_data[i, ]
  if (period$start * 100 <= max_time) {
    mid_x <- (period$end * 100 + min(period$start * 100, max_time)) / 2
    p <- p +
      annotate("text",
               x = mid_x,
               y = length(tree_time$tip.label) + 5,
               label = period$abbri,
               size = 2.5, fontface = "bold",
               color = isme_colors$neutral)
  }
}

# 主题设置
p <- p +
  labs(x = "Time (Mya)") +
  theme(
    legend.position = c(0.15, 0.85),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6),
    axis.title.x = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(size = 6),
    plot.margin = margin(5, 5, 5, 5)
  )

## =========================================================
## 6. 保存
## =========================================================

cat("保存图片...\n")

ggsave("Figure1_ggtree_single.pdf", p,
       width = 85, height = 150, units = "mm", dpi = 600)

ggsave("Figure1_ggtree_single.png", p,
       width = 85, height = 150, units = "mm", dpi = 600)

ggsave("Figure1_ggtree_double.pdf", p,
       width = 180, height = 150, units = "mm", dpi = 600)

ggsave("Figure1_ggtree_double.png", p,
       width = 180, height = 150, units = "mm", dpi = 600)

cat("\n✓ 完成!\n")
cat("  - Figure1_ggtree_single.pdf/png\n")
cat("  - Figure1_ggtree_double.pdf/png\n")
