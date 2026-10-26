#!/usr/bin/env Rscript
# Release path adapter: data and outputs live in a separate user workspace.
.release_root <- Sys.getenv("THERMAL_PROJECT_ROOT")
if (!nzchar(.release_root)) stop("Set THERMAL_PROJECT_ROOT to the analysis workspace")
.release_root <- normalizePath(.release_root, winslash = "/", mustWork = TRUE)
.release_input <- function(name) {
  value <- Sys.getenv(name)
  if (!nzchar(value)) stop("Set required input environment variable: ", name)
  normalizePath(value, winslash = "/", mustWork = TRUE)
}
setwd(file.path(.release_root, "abundance_wca_np"))

## =========================================================
## calc_NP_WCA_ratio.R
## 计算全球海洋 NP/WCA 比例
##   方法A: 观测样本面积加权平均
##   方法B: RF模型预测全球网格 → 面积加权
## =========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ranger)
  library(ggplot2)
  library(scales)
  library(ncdf4)
})

## ---- 路径 ----
base_dir <- file.path(.release_root, "abundance_wca_np")
np_wca_file <- file.path(base_dir, "NP_WCA_final.xlsx")
geo_file    <- file.path(base_dir, "genome_lon_lat.xlsx")

biooracle_grid <- .release_input("THERMAL_BASELINE_GRID")
bathy_nc <- .release_input("THERMAL_BATHYMETRY_NC")

out_dir <- file.path(base_dir, "NP_WCA_global_ratio")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

nearshore_km <- 200
shelf_depth_m <- 200
grid_res_deg <- 0.25
nseed  <- 10
ntree  <- 1000
min_node <- 3

## ---- 辅助函数 ----
area_weight <- function(lat) cos(lat * pi / 180)

align_lonlat_to_grid <- function(df, res_deg,
                                 xlim = c(-180, 180),
                                 ylim = c(-90,  90)) {
  df %>%
    mutate(
      Longitude = floor((.data[["Longitude"]] - xlim[1]) / res_deg) * res_deg + xlim[1] + res_deg / 2,
      Latitude  = floor((.data[["Latitude"]]  - ylim[1]) / res_deg) * res_deg + ylim[1] + res_deg / 2
    )
}

## =========================================================
## 1) 读取数据并合并
## =========================================================
message("=== 读取 NP_WCA_final.xlsx 和 genome_lon_lat.xlsx ===")

np_wca <- read_excel(np_wca_file) %>%
  rename(sample = 1, NP_frac = 2, WCA_frac = 3) %>%
  mutate(
    NP_frac  = suppressWarnings(as.numeric(NP_frac)),
    WCA_frac = suppressWarnings(as.numeric(WCA_frac))
  )

geo <- read_excel(geo_file) %>%
  rename(sample = 1, Longitude = 2, Latitude = 3) %>%
  mutate(
    Longitude = suppressWarnings(as.numeric(Longitude)),
    Latitude  = suppressWarnings(as.numeric(Latitude))
  ) %>%
  filter(is.finite(Longitude), is.finite(Latitude))

np_wca <- np_wca %>% distinct(sample, .keep_all = TRUE)
geo    <- geo    %>% distinct(sample, .keep_all = TRUE)

dat <- inner_join(np_wca, geo, by = "sample") %>%
  filter(is.finite(NP_frac), is.finite(WCA_frac))

message("合并后样本数: ", nrow(dat))

# 检查 NP + WCA 是否 = 1，如果不是则重新归一化
dat <- dat %>%
  mutate(
    total = NP_frac + WCA_frac,
    NP_frac = NP_frac / total,
    WCA_frac = WCA_frac / total
  ) %>%
  select(-total)

not_sum_1 <- sum(abs((dat$NP_frac + dat$WCA_frac) - 1) > 0.001)
if (not_sum_1 > 0) {
  message("警告: ", not_sum_1, " 个样本的 NP+WCA != 1，已重新归一化")
}

message("经度范围: ", round(min(dat$Longitude), 2), " ~ ", round(max(dat$Longitude), 2))
message("纬度范围: ", round(min(dat$Latitude), 2), " ~ ", round(max(dat$Latitude), 2))

## =========================================================
## 方法A: 直接面积加权平均
## =========================================================
message("\n=== 方法A: 观测样本面积加权平均 ===")

dat_a <- dat %>%
  mutate(w = area_weight(Latitude))

np_weighted_mean  <- sum(dat_a$NP_frac  * dat_a$w) / sum(dat_a$w)
wca_weighted_mean <- sum(dat_a$WCA_frac * dat_a$w) / sum(dat_a$w)
np_wca_ratio_a    <- np_weighted_mean / wca_weighted_mean

cat("\n--- 方法A 结果 ---\n")
cat("样本数:          ", nrow(dat_a), "\n")
cat("NP_frac 加权均值: ", round(np_weighted_mean, 4), "\n")
cat("WCA_frac 加权均值:", round(wca_weighted_mean, 4), "\n")
cat("NP/WCA 比值:      ", round(np_wca_ratio_a, 2), "\n")
cat("NP 占比:          ", round(np_weighted_mean * 100, 2), "%\n")
cat("WCA 占比:         ", round(wca_weighted_mean * 100, 2), "%\n")

simple_mean_np  <- mean(dat$NP_frac)
simple_mean_wca <- mean(dat$WCA_frac)
cat("\n(参考: 简单平均 NP=", round(simple_mean_np, 4),
    ", WCA=", round(simple_mean_wca, 4),
    ", NP/WCA=", round(simple_mean_np / simple_mean_wca, 2), ")\n")

## 按纬度带统计
lat_band_a <- dat_a %>%
  mutate(lat_band = floor(Latitude / 10) * 10) %>%
  group_by(lat_band) %>%
  summarise(
    n = n(),
    NP_mean  = weighted.mean(NP_frac, w),
    WCA_mean = weighted.mean(WCA_frac, w),
    .groups = "drop"
  ) %>%
  mutate(NP_WCA_ratio = NP_mean / pmax(WCA_mean, 1e-10))

cat("\n纬度带统计 (方法A):\n")
print(lat_band_a, n = Inf)

## =========================================================
## 方法B: RF模型预测全球网格
## =========================================================
message("\n=== 方法B: RF模型预测全球网格 ===")

## 2) 为训练数据提取Bio-ORACLE环境变量
message("读取Bio-ORACLE基线网格...")
base_grid <- read_tsv(biooracle_grid, show_col_types = FALSE, progress = FALSE)
message("基线网格行数: ", nrow(base_grid))

base_grid <- base_grid %>%
  mutate(
    Longitude = suppressWarnings(as.numeric(Longitude)),
    Latitude  = suppressWarnings(as.numeric(Latitude))
  ) %>%
  filter(is.finite(Longitude), is.finite(Latitude))

pred_vars <- c("temp", "sal", "chl", "si", "no3", "po4", "o2")
missing_vars <- setdiff(pred_vars, names(base_grid))
if (length(missing_vars) > 0) {
  message("警告: 基线网格缺少变量: ", paste(missing_vars, collapse = ", "))
  pred_vars <- intersect(pred_vars, names(base_grid))
}
message("预测变量: ", paste(pred_vars, collapse = ", "))

## 提取bathymetry用于nearshore定义
message("提取bathymetry (depth)...")
nc <- nc_open(bathy_nc)
bathy_lon <- ncvar_get(nc, "lon")
bathy_lat <- ncvar_get(nc, "lat")
bathy_val <- ncvar_get(nc, "bathymetry_mean")
nc_close(nc)

bathy_df <- expand.grid(Longitude = bathy_lon, Latitude = bathy_lat) %>%
  mutate(depth_m = as.vector(bathy_val)) %>%
  filter(!is.na(depth_m))
message("Bathymetry网格点数: ", nrow(bathy_df))

## 3) 为每个样本点提取最近网格的环境变量
message("为样本点提取环境变量...")

base_grid_aligned <- base_grid %>%
  align_lonlat_to_grid(grid_res_deg) %>%
  group_by(Longitude, Latitude) %>%
  summarise(across(all_of(c(pred_vars, "dist_to_land_km")),
                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")

bathy_aligned <- bathy_df %>%
  align_lonlat_to_grid(grid_res_deg) %>%
  group_by(Longitude, Latitude) %>%
  summarise(depth_m = mean(depth_m, na.rm = TRUE), .groups = "drop")

base_grid_aligned <- base_grid_aligned %>%
  left_join(bathy_aligned, by = c("Longitude", "Latitude"))

depth_na_rate <- mean(is.na(base_grid_aligned$depth_m))
message("基线网格 depth_m NA率: ", round(depth_na_rate * 100, 1), "%")

dat_aligned <- dat %>%
  align_lonlat_to_grid(grid_res_deg)

dat_env <- dat_aligned %>%
  left_join(base_grid_aligned, by = c("Longitude", "Latitude"), suffix = c("", "_grid"))

na_rate <- mean(is.na(dat_env$temp))
message("环境变量NA率: ", round(na_rate * 100, 1), "%")

train_dat <- dat_env %>%
  select(sample, Longitude, Latitude, NP_frac, WCA_frac,
         all_of(pred_vars), dist_to_land_km, depth_m) %>%
  drop_na(NP_frac, WCA_frac, all_of(pred_vars))
message("有效训练样本: ", nrow(train_dat))

## ---- Shelf-associated samples: latitude & SST range ----
shelf_samples <- dat_env %>%
  filter(is.finite(dist_to_land_km), dist_to_land_km <= nearshore_km,
         is.finite(depth_m),         depth_m >= -shelf_depth_m,
         is.finite(temp),            is.finite(NP_frac))

cat("\n--- Shelf-associated NP samples (dist <= ", nearshore_km,
    " km & depth >= -", shelf_depth_m, " m) ---\n", sep = "")
cat("样本数:        ", nrow(shelf_samples), "\n")
cat("纬度范围:      ", round(min(shelf_samples$Latitude), 2), "° ~ ",
                       round(max(shelf_samples$Latitude), 2), "°\n", sep = "")
cat("SST 范围:      ", round(min(shelf_samples$temp), 2), " ~ ",
                       round(max(shelf_samples$temp), 2), " °C\n", sep = "")
cat("SST 中位/均值: ", round(median(shelf_samples$temp), 2), " / ",
                       round(mean(shelf_samples$temp), 2), " °C\n", sep = "")
cat("NP_frac 中位/均值: ", round(median(shelf_samples$NP_frac), 3), " / ",
                           round(mean(shelf_samples$NP_frac), 3), "\n", sep = "")

shelf_np_dom <- shelf_samples %>% filter(NP_frac > 0.5)
cat("  其中 NP_frac>0.5 (NP占优, n=", nrow(shelf_np_dom), "):\n", sep = "")
cat("    纬度: ", round(min(shelf_np_dom$Latitude), 2), "° ~ ",
                  round(max(shelf_np_dom$Latitude), 2), "°\n", sep = "")
cat("    SST:  ", round(min(shelf_np_dom$temp), 2), " ~ ",
                  round(max(shelf_np_dom$temp), 2), " °C\n", sep = "")

write_tsv(shelf_samples %>%
            select(sample, Longitude, Latitude, depth_m, dist_to_land_km,
                   temp, NP_frac, WCA_frac),
          file.path(out_dir, "shelf_associated_samples.tsv"))
message("Shelf samples saved: ", file.path(out_dir, "shelf_associated_samples.tsv"))

## 4) 训练RF模型 (预测NP_frac)
message("训练RF集成模型 (", nseed, " seeds)...")

train_df <- train_dat %>% select(y = NP_frac, all_of(pred_vars))

fit_one_rf <- function(train_df, mtry, ntree, min_node, seed) {
  set.seed(seed)
  ranger(
    y ~ .,
    data = train_df,
    num.trees = ntree,
    mtry = mtry,
    min.node.size = min_node,
    splitrule = "variance"
  )
}

predict_ensemble <- function(train_df, newX, nseed = 10, ntree = 1000, min_node = 3) {
  p <- ncol(newX)
  mtry0 <- max(1, floor(sqrt(p)))
  mtry0 <- min(mtry0, p)
  pred_mat <- matrix(NA_real_, nrow = nrow(newX), ncol = nseed)
  for (i in seq_len(nseed)) {
    fit <- fit_one_rf(train_df, mtry = mtry0, ntree = ntree, min_node = min_node, seed = i)
    pred_mat[, i] <- predict(fit, data = newX)$predictions
  }
  rowMeans(pred_mat, na.rm = TRUE)
}

## OOB性能评估
set.seed(42)
rf_oob <- ranger(
  y ~ .,
  data = train_df,
  num.trees = ntree,
  min.node.size = min_node,
  splitrule = "variance",
  importance = "permutation"
)
cat("\nRF OOB R²: ", round(rf_oob$r.squared, 4), "\n")
cat("RF OOB RMSE: ", round(sqrt(rf_oob$prediction.error), 4), "\n")

## 变量重要性分析
var_imp <- importance(rf_oob)
var_imp_sorted <- sort(var_imp, decreasing = TRUE)
cat("\n--- 变量重要性 (permutation importance) ---\n")
for (v in names(var_imp_sorted)) {
  cat(sprintf("  %-6s: %.6f\n", v, var_imp_sorted[v]))
}
cat("最重要的预测因子: ", names(var_imp_sorted)[1], "\n")

imp_df <- tibble(
  variable = factor(names(var_imp_sorted), levels = rev(names(var_imp_sorted))),
  importance = as.numeric(var_imp_sorted)
)

## ---- 方向4: 近岸 vs 远洋分别的变量重要性 ----
## 拆分训练集，分别再训一个 RF，提取 permutation importance
fit_imp_subset <- function(sub_df, label, min_n = 15) {
  sub_train <- sub_df %>%
    select(y = NP_frac, all_of(pred_vars))
  n <- nrow(sub_train)
  cat("\n--- ", label, " 子集 RF (n=", n, ") ---\n", sep = "")
  if (n < min_n) {
    message("  样本量 < ", min_n, "，跳过 ", label, " 子集 RF")
    return(NULL)
  }
  set.seed(42)
  fit <- ranger(
    y ~ .,
    data = sub_train,
    num.trees = ntree,
    min.node.size = min_node,
    splitrule = "variance",
    importance = "permutation"
  )
  cat("  OOB R²:   ", round(fit$r.squared, 4), "\n")
  cat("  OOB RMSE: ", round(sqrt(fit$prediction.error), 4), "\n")
  imp_sorted <- sort(importance(fit), decreasing = TRUE)
  for (v in names(imp_sorted)) {
    cat(sprintf("  %-6s: %.6f\n", v, imp_sorted[v]))
  }
  tibble(
    subset = label,
    n = n,
    OOB_R2 = round(fit$r.squared, 4),
    variable = names(imp_sorted),
    importance = as.numeric(imp_sorted)
  )
}

nearshore_train <- train_dat %>%
  filter(is.finite(dist_to_land_km), dist_to_land_km <= nearshore_km,
         is.finite(depth_m),         depth_m >= -shelf_depth_m)

offshore_train <- train_dat %>%
  filter(!(is.finite(dist_to_land_km) & dist_to_land_km <= nearshore_km &
           is.finite(depth_m)         & depth_m >= -shelf_depth_m))

imp_nearshore <- fit_imp_subset(nearshore_train, "nearshore")
imp_offshore  <- fit_imp_subset(offshore_train,  "offshore")

imp_subset_df <- bind_rows(
  imp_df %>% mutate(subset = "global", n = nrow(train_dat)) %>%
    select(subset, n, variable, importance) %>%
    mutate(variable = as.character(variable)),
  imp_nearshore %>% select(subset, n, variable, importance),
  imp_offshore  %>% select(subset, n, variable, importance)
)

write_tsv(imp_subset_df,
          file.path(out_dir, "variable_importance_by_subset.tsv"))
message("Subset variable importance saved: ",
        file.path(out_dir, "variable_importance_by_subset.tsv"))

## 5) 全球海洋网格预测
message("预测全球海洋网格...")

ocean_grid <- base_grid_aligned %>%
  filter(complete.cases(across(all_of(pred_vars))))
message("全球海洋网格点数: ", nrow(ocean_grid))

ocean_newX <- ocean_grid %>% select(all_of(pred_vars))
ocean_grid$pred_NP_frac <- predict_ensemble(train_df, ocean_newX,
                                             nseed = nseed, ntree = ntree, min_node = min_node)

# 限制在 [0, 1]
ocean_grid$pred_NP_frac <- pmin(pmax(ocean_grid$pred_NP_frac, 0), 1)
ocean_grid$pred_WCA_frac <- 1 - ocean_grid$pred_NP_frac

## 6) 面积加权全球比例
ocean_grid <- ocean_grid %>%
  mutate(w = area_weight(Latitude))

np_global_b  <- sum(ocean_grid$pred_NP_frac  * ocean_grid$w) / sum(ocean_grid$w)
wca_global_b <- sum(ocean_grid$pred_WCA_frac * ocean_grid$w) / sum(ocean_grid$w)
np_wca_ratio_b <- np_global_b / wca_global_b

cat("\n--- 方法B 结果 (全球海洋) ---\n")
cat("网格点数:          ", nrow(ocean_grid), "\n")
cat("NP_frac 加权均值:  ", round(np_global_b, 4), "\n")
cat("WCA_frac 加权均值: ", round(wca_global_b, 4), "\n")
cat("NP/WCA 比值:       ", round(np_wca_ratio_b, 2), "\n")
cat("NP 占比:           ", round(np_global_b * 100, 2), "%\n")
cat("WCA 占比:          ", round(wca_global_b * 100, 2), "%\n")

## 近岸 vs 远洋 (nearshore = dist <= 200km AND depth >= -200m)
if ("dist_to_land_km" %in% names(ocean_grid) && "depth_m" %in% names(ocean_grid)) {
  nearshore <- ocean_grid %>%
    filter(dist_to_land_km <= nearshore_km, depth_m >= -shelf_depth_m)
  offshore  <- ocean_grid %>%
    filter(!(dist_to_land_km <= nearshore_km & depth_m >= -shelf_depth_m))

  np_nearshore  <- sum(nearshore$pred_NP_frac * nearshore$w) / sum(nearshore$w)
  wca_nearshore <- sum(nearshore$pred_WCA_frac * nearshore$w) / sum(nearshore$w)

  np_offshore  <- sum(offshore$pred_NP_frac * offshore$w) / sum(offshore$w)
  wca_offshore <- sum(offshore$pred_WCA_frac * offshore$w) / sum(offshore$w)

  cat("\n--- 近岸 (<=", nearshore_km, "km & depth >= -", shelf_depth_m, "m) ---\n")
  cat("网格点数: ", nrow(nearshore), "\n")
  cat("NP: ", round(np_nearshore * 100, 2), "%, WCA: ", round(wca_nearshore * 100, 2), "%\n")
  cat("NP/WCA: ", round(np_nearshore / wca_nearshore, 2), "\n")

  cat("\n--- 远洋 (其余) ---\n")
  cat("网格点数: ", nrow(offshore), "\n")
  cat("NP: ", round(np_offshore * 100, 2), "%, WCA: ", round(wca_offshore * 100, 2), "%\n")
  cat("NP/WCA: ", round(np_offshore / wca_offshore, 2), "\n")
}

## 纬度带统计
lat_band_b <- ocean_grid %>%
  mutate(lat_band = floor(Latitude / 10) * 10) %>%
  group_by(lat_band) %>%
  summarise(
    n = n(),
    NP_mean  = weighted.mean(pred_NP_frac, w),
    WCA_mean = weighted.mean(pred_WCA_frac, w),
    .groups = "drop"
  ) %>%
  mutate(NP_WCA_ratio = NP_mean / pmax(WCA_mean, 1e-10))

cat("\n纬度带统计 (方法B):\n")
print(lat_band_b, n = Inf)

## 把关键统计量写到一个 summary tsv
key_stats <- data.frame(
  metric = c("n_train", "OOB_R2", "OOB_RMSE",
             "most_important_var",
             "NP_frac_nearshore_pct", "NP_frac_offshore_pct",
             "NP_frac_global_pct"),
  value = c(nrow(train_dat),
            round(rf_oob$r.squared, 4),
            round(sqrt(rf_oob$prediction.error), 4),
            names(var_imp_sorted)[1],
            round(np_nearshore * 100, 2),
            round(np_offshore * 100, 2),
            round(np_global_b * 100, 2))
)
write_tsv(key_stats, file.path(out_dir, "key_stats_for_manuscript.tsv"))
message("Key stats saved: ", file.path(out_dir, "key_stats_for_manuscript.tsv"))

## 保存变量重要性
write_tsv(imp_df %>% arrange(desc(importance)),
          file.path(out_dir, "variable_importance.tsv"))
message("Variable importance saved: ", file.path(out_dir, "variable_importance.tsv"))

## =========================================================
## 7) 两种方法对比总结
## =========================================================
cat("\n\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("  全球 NP vs WCA 比例总结\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("\n方法A (观测样本面积加权):\n")
cat("  NP : WCA = ", round(np_weighted_mean * 100, 1), "% : ",
    round(wca_weighted_mean * 100, 1), "%\n")
cat("  NP/WCA = ", round(np_wca_ratio_a, 1), "\n")

cat("\n方法B (RF全球网格预测):\n")
cat("  NP : WCA = ", round(np_global_b * 100, 1), "% : ",
    round(wca_global_b * 100, 1), "%\n")
cat("  NP/WCA = ", round(np_wca_ratio_b, 1), "\n")

## =========================================================
## 8) 可视化
## =========================================================
message("\n=== 生成图表 ===")

suppressPackageStartupMessages(library(patchwork))

## Fig 1: 全球地图 — RF预测的NP_frac
world_df <- ggplot2::map_data("world")

p_map <- ggplot() +
  geom_tile(data = ocean_grid,
            aes(x = Longitude, y = Latitude, fill = pred_NP_frac),
            width = grid_res_deg, height = grid_res_deg) +
  geom_polygon(data = world_df, aes(x = long, y = lat, group = group),
               fill = "grey92", color = "grey55", linewidth = 0.15) +
  scale_fill_gradient2(low = "#58D68D", mid = "white", high = "#E8A07A",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "NP fraction") +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-70, 80), expand = FALSE) +
  labs(title = "Global NP vs WCA relative abundance (RF prediction)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "#D6EAF8", color = NA),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4),
    plot.title = element_text(face = "bold", size = 11)
  )

## Fig 2: 观测点叠加
p_obs <- ggplot() +
  geom_tile(data = ocean_grid,
            aes(x = Longitude, y = Latitude, fill = pred_NP_frac),
            width = grid_res_deg, height = grid_res_deg, alpha = 0.6) +
  geom_polygon(data = world_df, aes(x = long, y = lat, group = group),
               fill = "grey92", color = "grey55", linewidth = 0.15) +
  geom_point(data = dat, aes(x = Longitude, y = Latitude, color = NP_frac),
             size = 1.5, shape = 16) +
  scale_fill_gradient2(low = "#58D68D", mid = "white", high = "#E8A07A",
                       midpoint = 0.5, limits = c(0, 1),
                       name = "NP frac\n(grid)") +
  scale_color_gradient2(low = "#58D68D", mid = "white", high = "#E8A07A",
                        midpoint = 0.5, limits = c(0, 1),
                        name = "NP frac\n(obs)") +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-70, 80), expand = FALSE) +
  labs(title = "Observed samples overlaid on RF prediction",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.background = element_rect(fill = "#D6EAF8", color = NA),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4),
    plot.title = element_text(face = "bold", size = 11)
  )

## Fig 3: 纬度梯度对比
lat_compare <- bind_rows(
  lat_band_a %>% mutate(method = "A: Observed"),
  lat_band_b %>% mutate(method = "B: RF predicted")
)

p_lat <- ggplot(lat_compare, aes(x = lat_band, y = WCA_mean * 100,
                                  color = method, fill = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("A: Observed" = "#2E86C1", "B: RF predicted" = "#E74C3C")) +
  labs(x = "Latitude", y = "WCA fraction (%)",
       title = "WCA fraction by latitude band",
       color = NULL, fill = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")

## Fig 4: 变量重要性
p_imp <- ggplot(imp_df, aes(x = variable, y = importance)) +
  geom_col(fill = "#2E86C1", width = 0.6) +
  coord_flip() +
  labs(x = NULL, y = "Permutation importance",
       title = "Variable importance for NP fraction prediction") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

## Fig 4b: 近岸 vs 远洋变量重要性对比
if (!is.null(imp_nearshore) || !is.null(imp_offshore)) {
  imp_subset_plot <- imp_subset_df %>%
    mutate(subset = factor(subset, levels = c("global", "nearshore", "offshore")),
           variable = factor(variable, levels = rev(names(var_imp_sorted))))

  p_imp_subset <- ggplot(imp_subset_plot,
                          aes(x = variable, y = importance, fill = subset)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    coord_flip() +
    scale_fill_manual(values = c(global = "#5D6D7E",
                                  nearshore = "#E67E22",
                                  offshore  = "#2E86C1")) +
    labs(x = NULL, y = "Permutation importance",
         fill = NULL,
         title = "Variable importance: nearshore vs offshore") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")

  out_imp_pdf <- file.path(out_dir, "variable_importance_by_subset.pdf")
  ggsave(out_imp_pdf, p_imp_subset, width = 8, height = 5, dpi = 300, bg = "white")
  message("Subset importance PDF saved: ", out_imp_pdf)
}

## Fig 5: 方法对比柱状图
compare_df <- tibble(
  method = c(paste0("A: Observed\n(n=", nrow(dat_a), ")"), paste0("B: RF predicted\n(n=", nrow(ocean_grid), ")")),
  NP = c(np_weighted_mean, np_global_b),
  WCA = c(wca_weighted_mean, wca_global_b)
) %>%
  pivot_longer(cols = c(NP, WCA), names_to = "group", values_to = "fraction")

p_bar <- ggplot(compare_df, aes(x = method, y = fraction * 100, fill = group)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = c(NP = "#E8A07A", WCA = "#58D68D")) +
  labs(x = NULL, y = "Fraction (%)", fill = NULL,
       title = "NP vs WCA: two methods comparison") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

## 组合输出
combined <- (p_map / p_obs) | (p_lat / p_imp / p_bar)
combined <- combined +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 13)))

out_pdf <- file.path(out_dir, "NP_WCA_global_ratio_comparison_v2.pdf")
ggsave(out_pdf, combined, width = 16, height = 12, dpi = 300, bg = "white")
message("PDF saved: ", out_pdf)

## 单独保存全球地图
out_map_pdf <- file.path(out_dir, "NP_WCA_global_map_RF_v2.pdf")
if (file.exists(out_map_pdf)) {
  message("Map PDF already exists, skipping: ", out_map_pdf)
} else {
  ggsave(out_map_pdf, p_map, width = 12, height = 6, dpi = 300, bg = "white")
  message("Map PDF saved: ", out_map_pdf)
}

## 保存结果表
write_tsv(ocean_grid %>% select(Longitude, Latitude, pred_NP_frac, pred_WCA_frac),
          file.path(out_dir, "global_ocean_NP_WCA_predicted.tsv"))

write_tsv(lat_compare, file.path(out_dir, "latband_NP_WCA_comparison.tsv"))

## =========================================================
## 9) SSP585 未来预测
## =========================================================
message("\n=== SSP585 未来预测 (2080-2100) ===")

future_grid_file <- .release_input("THERMAL_FUTURE_GRID")

if (!file.exists(future_grid_file)) {
  message("未找到SSP585未来网格文件，跳过: ", future_grid_file)
} else {
  fut_grid <- read_tsv(future_grid_file, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      Longitude = suppressWarnings(as.numeric(Longitude)),
      Latitude  = suppressWarnings(as.numeric(Latitude))
    ) %>%
    filter(is.finite(Longitude), is.finite(Latitude))
  message("未来网格行数: ", nrow(fut_grid))

  fut_aligned <- fut_grid %>%
    align_lonlat_to_grid(grid_res_deg) %>%
    group_by(Longitude, Latitude) %>%
    summarise(across(all_of(intersect(pred_vars, names(fut_grid))),
                     ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  fut_ocean <- fut_aligned %>%
    filter(complete.cases(across(all_of(pred_vars))))
  message("未来有效网格点: ", nrow(fut_ocean))

  fut_newX <- fut_ocean %>% select(all_of(pred_vars))
  fut_ocean$pred_NP_frac <- predict_ensemble(train_df, fut_newX,
                                              nseed = nseed, ntree = ntree, min_node = min_node)
  fut_ocean$pred_NP_frac <- pmin(pmax(fut_ocean$pred_NP_frac, 0), 1)
  fut_ocean$pred_WCA_frac <- 1 - fut_ocean$pred_NP_frac
  fut_ocean <- fut_ocean %>% mutate(w = area_weight(Latitude))

  np_fut  <- sum(fut_ocean$pred_NP_frac  * fut_ocean$w) / sum(fut_ocean$w)
  wca_fut <- sum(fut_ocean$pred_WCA_frac * fut_ocean$w) / sum(fut_ocean$w)

  cat("\n--- SSP585 2080-2100 全球预测 ---\n")
  cat("网格点数:          ", nrow(fut_ocean), "\n")
  cat("NP_frac 加权均值:  ", round(np_fut, 4), "\n")
  cat("WCA_frac 加权均值: ", round(wca_fut, 4), "\n")
  cat("NP/WCA 比值:       ", round(np_fut / wca_fut, 2), "\n")
  cat("NP 占比:           ", round(np_fut * 100, 2), "%\n")
  cat("WCA 占比:          ", round(wca_fut * 100, 2), "%\n")

  ## 加入dist_to_land_km和depth_m做近岸/远洋分析
  fut_ocean <- fut_ocean %>%
    left_join(base_grid_aligned %>% select(Longitude, Latitude, dist_to_land_km, depth_m),
              by = c("Longitude", "Latitude"))

  if ("dist_to_land_km" %in% names(fut_ocean) && "depth_m" %in% names(fut_ocean)) {
    ns_fut <- fut_ocean %>%
      filter(is.finite(dist_to_land_km), dist_to_land_km <= nearshore_km,
             is.finite(depth_m), depth_m >= -shelf_depth_m)
    os_fut <- fut_ocean %>%
      filter(!(is.finite(dist_to_land_km) & dist_to_land_km <= nearshore_km &
               is.finite(depth_m) & depth_m >= -shelf_depth_m))

    np_ns_fut  <- sum(ns_fut$pred_NP_frac * ns_fut$w) / sum(ns_fut$w)
    wca_ns_fut <- sum(ns_fut$pred_WCA_frac * ns_fut$w) / sum(ns_fut$w)
    np_os_fut  <- sum(os_fut$pred_NP_frac * os_fut$w) / sum(os_fut$w)
    wca_os_fut <- sum(os_fut$pred_WCA_frac * os_fut$w) / sum(os_fut$w)

    cat("\n--- SSP585 近岸 (<=", nearshore_km, "km & depth >= -", shelf_depth_m, "m) ---\n")
    cat("NP: ", round(np_ns_fut * 100, 2), "%, WCA: ", round(wca_ns_fut * 100, 2), "%\n")
    cat("NP/WCA: ", round(np_ns_fut / wca_ns_fut, 2), "\n")

    cat("\n--- SSP585 远洋 (其余) ---\n")
    cat("NP: ", round(np_os_fut * 100, 2), "%, WCA: ", round(wca_os_fut * 100, 2), "%\n")
    cat("NP/WCA: ", round(np_os_fut / wca_os_fut, 2), "\n")
  }

  ## 纬度带
  lat_band_fut <- fut_ocean %>%
    mutate(lat_band = floor(Latitude / 10) * 10) %>%
    group_by(lat_band) %>%
    summarise(
      n = n(),
      NP_mean  = weighted.mean(pred_NP_frac, w),
      WCA_mean = weighted.mean(pred_WCA_frac, w),
      .groups = "drop"
    ) %>%
    mutate(NP_WCA_ratio = NP_mean / pmax(WCA_mean, 1e-10))

  cat("\n纬度带统计 (SSP585 2080-2100):\n")
  print(lat_band_fut, n = Inf)

  ## Delta: 未来 - 当前
  delta_df <- fut_ocean %>%
    select(Longitude, Latitude, pred_NP_frac_fut = pred_NP_frac) %>%
    inner_join(ocean_grid %>% select(Longitude, Latitude, pred_NP_frac_cur = pred_NP_frac),
               by = c("Longitude", "Latitude")) %>%
    mutate(delta_NP = pred_NP_frac_fut - pred_NP_frac_cur)

  cat("\n--- Delta (SSP585 - Current) ---\n")
  cat("Delta NP_frac 均值: ", round(mean(delta_df$delta_NP, na.rm = TRUE), 4), "\n")
  cat("Delta NP_frac 中位数: ", round(median(delta_df$delta_NP, na.rm = TRUE), 4), "\n")
  cat("NP增加的网格占比: ", round(mean(delta_df$delta_NP > 0, na.rm = TRUE) * 100, 1), "%\n")
  cat("WCA增加的网格占比: ", round(mean(delta_df$delta_NP < 0, na.rm = TRUE) * 100, 1), "%\n")

  ## 总结对比
  cat("\n\n")
  cat("=" |> rep(60) |> paste(collapse = ""), "\n")
  cat("  Current vs SSP585 全球 NP/WCA 对比\n")
  cat("=" |> rep(60) |> paste(collapse = ""), "\n")
  cat("\n当前 (baseline):  NP ", round(np_global_b * 100, 1), "% : WCA ",
      round(wca_global_b * 100, 1), "%  (NP/WCA = ", round(np_global_b / wca_global_b, 2), ")\n")
  cat("SSP585 2080-2100: NP ", round(np_fut * 100, 1), "% : WCA ",
      round(wca_fut * 100, 1), "%  (NP/WCA = ", round(np_fut / wca_fut, 2), ")\n")
  cat("变化: NP ", ifelse(np_fut > np_global_b, "+", ""),
      round((np_fut - np_global_b) * 100, 2), " pp\n")

  ## =========================================================
  ## 10) SSP585 可视化
  ## =========================================================
  ## 未来地图
  p_map_fut <- ggplot() +
    geom_tile(data = fut_ocean,
              aes(x = Longitude, y = Latitude, fill = pred_NP_frac),
              width = grid_res_deg, height = grid_res_deg) +
    geom_polygon(data = world_df, aes(x = long, y = lat, group = group),
                 fill = "grey92", color = "grey55", linewidth = 0.15) +
    scale_fill_gradient2(low = "#58D68D", mid = "white", high = "#E8A07A",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "NP fraction") +
    coord_quickmap(xlim = c(-180, 180), ylim = c(-70, 80), expand = FALSE) +
    labs(title = "SSP585 2080-2100: NP vs WCA (RF prediction)",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "#D6EAF8", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4),
      plot.title = element_text(face = "bold", size = 11)
    )

  ## Delta地图
  q98 <- quantile(abs(delta_df$delta_NP), 0.98, na.rm = TRUE)

  p_delta <- ggplot() +
    geom_tile(data = delta_df,
              aes(x = Longitude, y = Latitude, fill = delta_NP),
              width = grid_res_deg, height = grid_res_deg) +
    geom_polygon(data = world_df, aes(x = long, y = lat, group = group),
                 fill = "grey92", color = "grey55", linewidth = 0.15) +
    scale_fill_gradient2(low = "#58D68D", mid = "white", high = "#E8A07A",
                         midpoint = 0, limits = c(-q98, q98), oob = squish,
                         name = expression(Delta * " NP frac")) +
    coord_quickmap(xlim = c(-180, 180), ylim = c(-70, 80), expand = FALSE) +
    labs(title = "Change in NP fraction (SSP585 2080-2100 minus Current)",
         subtitle = paste0("Green = WCA increases, Orange = NP increases"),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "#D6EAF8", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4),
      plot.title = element_text(face = "bold", size = 11)
    )

  ## 纬度带对比
  lat_compare_fut <- bind_rows(
    lat_band_b %>% mutate(period = "Current"),
    lat_band_fut %>% mutate(period = "SSP585 2080-2100")
  )

  p_lat_fut <- ggplot(lat_compare_fut, aes(x = lat_band, y = NP_mean * 100,
                                             color = period)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Current" = "#2E86C1", "SSP585 2080-2100" = "#E74C3C")) +
    labs(x = "Latitude", y = "NP fraction (%)",
         title = "NP fraction by latitude: Current vs SSP585",
         color = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")

  ## 组合
  combined_ssp <- (p_map / p_map_fut) | (p_delta / p_lat_fut)
  combined_ssp <- combined_ssp +
    plot_annotation(
      title = "NP vs WCA: Current vs SSP585 2080-2100",
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.tag = element_text(face = "bold", size = 13)
      )
    )

  out_ssp_pdf <- file.path(out_dir, "NP_WCA_current_vs_SSP585.pdf")
  ggsave(out_ssp_pdf, combined_ssp, width = 16, height = 12, dpi = 300, bg = "white")
  message("SSP585 PDF saved: ", out_ssp_pdf)

  ## 保存数据
  write_tsv(fut_ocean %>% select(Longitude, Latitude, pred_NP_frac, pred_WCA_frac),
            file.path(out_dir, "global_ocean_NP_WCA_SSP585_predicted.tsv"))
  write_tsv(delta_df, file.path(out_dir, "global_ocean_NP_WCA_delta_SSP585.tsv"))
  write_tsv(lat_compare_fut, file.path(out_dir, "latband_NP_WCA_current_vs_SSP585.tsv"))
}

message("\nALL DONE")
