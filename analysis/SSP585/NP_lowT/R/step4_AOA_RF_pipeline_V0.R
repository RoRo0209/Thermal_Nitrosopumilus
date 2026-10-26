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
setwd(thermal_path("SSP585/NP_lowT"))

## =========================================================
## AOA macroecology RF pipeline (CLEAN & FIXED)
## VIF -> (optional) RFE -> tuning -> ensemble -> importance -> PDP
## Input : AOA_model_table_V0_final.tsv
## Output:
##   AOA_model_table_V0_VIF.tsv
##   RFE_result.rds (optional) + RFE_selected_features.txt
##   RF_tuned_model.rds + RF_grid_results.csv
##   RF_importance.tsv + RF_importance.png
##   RF_ensemble_predictions.tsv
##   PDP_*.png
## =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(car)
  library(caret)
  library(ranger)
  library(tibble)
  library(pdp)
  library(ggplot2)
})

## -----------------------------
## 0) Paths (EDIT)
## -----------------------------
setwd(thermal_path("SSP585/NP_lowT/AOA_RF_V0_outputs"))
in_file_v0 <- "AOA_model_table_V0_final.tsv"

## -----------------------------
## 1) Read data + keep numeric
## -----------------------------
df0 <- read_tsv(in_file_v0, show_col_types = FALSE)

# 强制转 numeric（避免字符列混入）
df0 <- df0 %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

stopifnot("y" %in% names(df0))

pred_cols0 <- setdiff(colnames(df0), c("y","Longitude","Latitude"))

message("Input rows: ", nrow(df0), " ; predictors: ", length(pred_cols0))

## =========================================================
## Step 1: VIF filtering  (complete cases only)
## =========================================================
message("\n## Step1: VIF filtering")

df_vif <- df0 %>%
  select(y, all_of(pred_cols0), any_of(c("Longitude","Latitude"))) %>%
  drop_na(y, all_of(pred_cols0))

message("Rows for VIF (complete cases): ", nrow(df_vif))
if (nrow(df_vif) < 20) warning("Too few rows for reliable VIF (complete cases < 20).")

vif_loop <- function(df, thresh = 10) {
  vars <- setdiff(colnames(df), c("y","Longitude","Latitude"))
  repeat {
    if (length(vars) < 2) break
    fit <- lm(reformulate(vars, response = "y"), data = df)
    v <- car::vif(fit)
    if (length(v) < 2) break
    if (max(v, na.rm = TRUE) < thresh) break
    drop_var <- names(which.max(v))
    message("Dropping due to high VIF: ", drop_var, " (", round(max(v), 2), ")")
    df <- df %>% select(-all_of(drop_var))
    vars <- setdiff(vars, drop_var)
  }
  df
}

df_vif2 <- vif_loop(df_vif, thresh = 10)
write_tsv(df_vif2, "AOA_model_table_V0_VIF.tsv")
message("Wrote: AOA_model_table_V0_VIF.tsv (rows=", nrow(df_vif2), ", cols=", ncol(df_vif2), ")")

pred_cols1 <- setdiff(colnames(df_vif2), c("y","Longitude","Latitude"))

## =========================================================
## Step 2: RFE feature selection (optional; only if enough rows)
## =========================================================
message("\n## Step2: RFE feature selection")

# RFE 建议也用 complete cases
df_rfe <- df_vif2 %>% select(y, all_of(pred_cols1)) %>% drop_na()

message("Rows for RFE (complete cases): ", nrow(df_rfe))

if (nrow(df_rfe) < 200) {  # CHANGED: skip RFE for small datasets
  message("Rows < 200 -> SKIP RFE, use all predictors after VIF.")
  best_vars <- pred_cols1
  writeLines(best_vars, "RFE_selected_features.txt")
} else {
  x <- df_rfe[, pred_cols1]
  y <- df_rfe$y
  
  ctrl <- rfeControl(
    functions = caretFuncs,
    method = "cv",
    number = 10,
    verbose = TRUE
  )
  
  set.seed(1)
  rfe_fit <- rfe(
    x = x, y = y,
    sizes = seq(2, length(pred_cols1)),
    rfeControl = ctrl,
    method = "ranger",
    num.trees = 500
  )
  
  saveRDS(rfe_fit, "RFE_result.rds")
  best_vars <- predictors(rfe_fit)
  writeLines(best_vars, "RFE_selected_features.txt")
  
  message("Selected features (RFE): ", paste(best_vars, collapse = ", "))
}

## =========================================================
## Step 3: Hyperparameter tuning (caret + ranger) on best_vars
## =========================================================
message("\n## Step3: RF tuning")

best_vars <- readLines("RFE_selected_features.txt")

df_rf <- df_vif2 %>%
  select(y, all_of(best_vars)) %>%
  drop_na()

message("Rows for tuning (complete cases): ", nrow(df_rf))
if (nrow(df_rf) < 30) stop("Too few rows for tuning/training (<30). Fix NA first.")

ctrl <- trainControl(method = "cv", number = 10)

p <- length(best_vars)
mtry_grid <- unique(pmax(1, round(c(sqrt(p), p/3, p/2, 2*p/3, p))))
grid <- expand.grid(
  mtry = mtry_grid,
  splitrule = "variance",
  min.node.size = c(1, 3, 5, 10)
)

set.seed(1)
rf_tuned <- train(
  y ~ .,
  data = df_rf,
  method = "ranger",
  trControl = ctrl,
  tuneGrid = grid,
  num.trees = 1000,
  importance = "permutation"
)

saveRDS(rf_tuned, "RF_tuned_model.rds")
write.csv(rf_tuned$results, "RF_grid_results.csv", row.names = FALSE)

message("Best tune:")
print(rf_tuned$bestTune)

## =========================================================
## Step 3.5: Importance (AFTER rf_tuned exists)
## =========================================================
message("\n## Step3.5: variable importance")

vi <- varImp(rf_tuned, scale = TRUE)
vi_df <- vi$importance %>%
  tibble::rownames_to_column("var") %>%
  arrange(desc(Overall))

write_tsv(vi_df, "RF_importance.tsv")
message("Wrote: RF_importance.tsv")

# 保存 importance 图
png("RF_importance.png", width = 1600, height = 1200, res = 200)
plot(vi, top = min(20, nrow(vi$importance)))
dev.off()
message("Wrote: RF_importance.png")

## =========================================================
## Step 4: 10-seed ensemble (train on full df_rf; predict on same df_rf)
## =========================================================
message("\n## Step4: 10-seed ensemble")

pred_mat <- matrix(NA_real_, nrow = nrow(df_rf), ncol = 10)

for (i in 1:10) {
  set.seed(i)
  fit <- ranger(
    y ~ .,
    data = df_rf,
    num.trees = 1000,
    mtry = rf_tuned$bestTune$mtry,
    min.node.size = rf_tuned$bestTune$min.node.size,
    importance = "permutation"
  )
  pred_mat[, i] <- predict(fit, data = df_rf)$predictions
}

pred_mean <- rowMeans(pred_mat)
pred_sd   <- apply(pred_mat, 1, sd)
pred_cv_abs <- pred_sd / pmax(abs(pred_mean), 1e-8)

out_pred <- tibble(
  pred_mean = pred_mean,
  pred_sd   = pred_sd,
  pred_cv_abs = pred_cv_abs
)

write_tsv(out_pred, "RF_ensemble_predictions.tsv")
message("Wrote: RF_ensemble_predictions.tsv")

## =========================================================
## Step 5: PDP plots (for vars with Overall>0)
## =========================================================
message("\n## Step5: PDP plots")

vars_pdp <- vi_df %>% filter(Overall > 0) %>% pull(var)
x_train <- df_rf %>% select(-y)

for (v in vars_pdp) {
  p_obj <- partial(
    object = rf_tuned$finalModel,
    pred.var = v,
    train = x_train,
    grid.resolution = 50
  )
  
  g <- autoplot(p_obj) +
    theme_bw() +
    labs(title = paste0("PDP: ", v), y = "Predicted y")
  
  ggsave(filename = paste0("PDP_", v, ".png"), plot = g, width = 6, height = 4, dpi = 200)
}

message("\nALL DONE ✅")
