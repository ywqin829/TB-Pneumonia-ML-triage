## Prediction utilities for the local research calculator. No public deployment.
load_tb_bundle <- function(root) {
  suppressPackageStartupMessages({
    library(caret)
    library(recipes)
    library(themis)
    library(xgboost)
  })
  source(file.path(root, "cache_guard.R"), local = TRUE)
  root <- rerun_resolve_root(root)
  rerun_verify_manifest(root)
  lock <- readRDS(file.path(root, "results/analysis_lock.rds"))
  stopifnot(identical(lock$winner, "xgbTree"),
            setequal(lock$panel, c("MCV", "PDW", "MONO.")))
  thresholds <- read.csv(file.path(root, "results/locked_thresholds.csv"),
                          stringsAsFactors = FALSE)
  low <- thresholds$threshold[thresholds$rule == "sensitivity_90"]
  high <- thresholds$threshold[thresholds$rule == "rulein_specificity_90"]
  balanced <- thresholds$threshold[thresholds$rule == "youden"]
  stopifnot(length(low) == 1, length(high) == 1, length(balanced) == 1,
            low > 0, high < 1, low < balanced, balanced < high,
            abs(low - 0.289648950099945) < 1e-12,
            abs(high - 0.63411009311676) < 1e-12,
            abs(balanced - 0.422575533390045) < 1e-12)
  prep <- readRDS(file.path(root, "results/prepared_data.rds"))
  ranges <- rbind(MCV = range(prep$train$MCV), PDW = range(prep$train$PDW),
                  MONO_percent = range(prep$train$MONO.) * 100)
  list(model = readRDS(file.path(root, "models/final_xgbTree.rds")),
       thresholds = c(balanced = balanced, low = low, high = high), ranges = ranges, root = root,
       features = lock$panel)
}

predict_tb_inputs <- function(bundle, mcv, pdw, mono_percent) {
  values <- c(mcv, pdw, mono_percent)
  if (length(values) != 3 || any(!is.finite(values))) stop("All three inputs must be finite numbers.")
  if (mcv <= 0 || pdw <= 0 || mono_percent < 0 || mono_percent > 100) {
    stop("MCV and PDW must be positive; MONO% must be between 0 and 100.")
  }
  input <- data.frame(MCV = mcv, PDW = pdw, check.names = FALSE)
  input[["MONO."]] <- mono_percent / 100
  score <- as.numeric(predict(bundle$model, newdata = input, type = "prob")$TB)
  stopifnot(length(score) == 1, is.finite(score), score >= 0, score <= 1)
  model_class <- if (score >= bundle$thresholds[["balanced"]]) "TB" else "Pneumonia"
  zone <- if (score < bundle$thresholds[["low"]]) "low" else {
    if (score < bundle$thresholds[["high"]]) "intermediate" else "high"
  }
  observed <- c(MCV = mcv, PDW = pdw, MONO_percent = mono_percent)
  outside <- names(observed)[observed < bundle$ranges[, 1] | observed > bundle$ranges[, 2]]
  list(score = score, model_class = model_class, zone = zone, outside_training_range = outside,
       raw_model_input = input)
}
