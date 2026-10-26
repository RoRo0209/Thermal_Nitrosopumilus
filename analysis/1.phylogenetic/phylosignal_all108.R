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
## 在全部108个AOA基因组上计算K和Pagel's lambda
## 目的：先验证OGT的系统发育信号，再用ASR定义生态型
## 解决循环论证问题
## =========================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

setwd(thermal_path("1.phylogenetic"))

## =========================================================
## 1. 读取数据
## =========================================================

# ML树（原始树，未经chronos转换）
tree_ml <- read.tree("gtdbtk.ar53.NP.treefile")
cat("原始ML树叶节点数:", length(tree_ml$tip.label), "\n")

# OGT数据（含全部108个基因组）
meta <- read.csv("ASR_clade_assignments_final.csv")
cat("元数据行数:", nrow(meta), "\n")
cat("OGT缺失值:", sum(is.na(meta$OGT)), "\n\n")

# 只保留树和元数据共有的基因组
common <- intersect(tree_ml$tip.label, meta$FileName)
tree_ml <- drop.tip(tree_ml, setdiff(tree_ml$tip.label, common))
meta <- meta[match(tree_ml$tip.label, meta$FileName), ]
cat("匹配后基因组数:", length(common), "\n\n")

## =========================================================
## 2. 准备三个分析子集
## =========================================================

# 子集A：全部108个基因组（无ASR筛选）
ogt_all <- setNames(meta$OGT, meta$FileName)
tree_all <- tree_ml

# 子集B：仅LowT + HighT clade（34个，ASR筛选后，原始分析）
clade_meta <- meta[meta$Clade %in% c("LowT", "HighT"), ]
ogt_clade <- setNames(clade_meta$OGT, clade_meta$FileName)
tree_clade <- keep.tip(tree_ml, names(ogt_clade))

# 子集C：全部108个基因组，使用chronos时间树
tree_chrono <- chronos(tree_ml, lambda = 1, quiet = TRUE)
ogt_all_chrono <- ogt_all  # OGT数据不变

cat("=== 分析子集 ===\n")
cat("子集A（全部，ML树）:", length(ogt_all), "个基因组\n")
cat("子集B（仅clade，ML树）:", length(ogt_clade), "个基因组\n")
cat("子集C（全部，时间树）:", length(ogt_all_chrono), "个基因组\n\n")

## =========================================================
## 3. 计算K和lambda
## =========================================================

run_phylosig <- function(tree, ogt, label) {
  cat("---", label, "---\n")

  K_res <- phylosig(tree, ogt, method = "K", test = TRUE, nsim = 1000)
  lam_res <- phylosig(tree, ogt, method = "lambda", test = TRUE)

  cat(sprintf("  Blomberg's K = %.4f  (p = %s)\n",
              K_res$K, format.pval(K_res$P, digits = 3)))
  cat(sprintf("  Pagel's λ   = %.4f  (p = %s)\n\n",
              lam_res$lambda, format.pval(lam_res$P, digits = 3)))

  list(K = K_res$K, K_p = K_res$P,
       lambda = lam_res$lambda, lambda_p = lam_res$P)
}

cat("=== 系统发育信号分析结果 ===\n\n")

res_A <- run_phylosig(tree_all,    ogt_all,        "A: 全部108个基因组（ML树）")
res_B <- run_phylosig(tree_clade,  ogt_clade,      "B: 仅LowT+HighT clade（ML树，原始分析）")
res_C <- run_phylosig(tree_chrono, ogt_all_chrono, "C: 全部108个基因组（时间树）")

## =========================================================
## 4. 汇总对比
## =========================================================

cat("=== 汇总对比 ===\n")
cat(sprintf("%-40s  %8s  %8s  %8s  %8s\n",
            "分析", "K", "K_p", "lambda", "lambda_p"))
cat(strrep("-", 80), "\n")
cat(sprintf("%-40s  %8.4f  %8s  %8.4f  %8s\n",
            "A: 全部108（ML树）",
            res_A$K, format.pval(res_A$K_p, digits=2),
            res_A$lambda, format.pval(res_A$lambda_p, digits=2)))
cat(sprintf("%-40s  %8.4f  %8s  %8.4f  %8s\n",
            "B: 仅clade 34（ML树，原始）",
            res_B$K, format.pval(res_B$K_p, digits=2),
            res_B$lambda, format.pval(res_B$lambda_p, digits=2)))
cat(sprintf("%-40s  %8.4f  %8s  %8.4f  %8s\n",
            "C: 全部108（时间树）",
            res_C$K, format.pval(res_C$K_p, digits=2),
            res_C$lambda, format.pval(res_C$lambda_p, digits=2)))

## =========================================================
