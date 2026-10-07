## Full revision rerun: same Boruta -> LASSO -> AUC -> combinations ->
## nine-algorithm route as the manuscript, with all learned steps inside
## training folds. Run after 00_prepare.R from the rerun directory.

suppressPackageStartupMessages({
  library(Boruta)
  library(glmnet)
  library(caret)
  library(recipes)
  library(themis)
  library(pROC)
})

set.seed(2026)
dir.create("results", showWarnings = FALSE)
dir.create("models", showWarnings = FALSE)
prep <- readRDS("results/prepared_data.rds")
train_data <- prep$train
candidate_features <- prep$candidate_features
positive <- "TB"
algorithms <- c("glm", "glmnet", "rpart", "ranger", "xgbTree",
                "svmRadial", "gbm", "nnet", "knn")
old_panel <- c("RBC", "PDW", "MONO.", "MPV", "LMR")
stopifnot(all(old_panel %in% names(train_data)))

auc_rank <- function(y, score) {
  pos <- y == positive
  npos <- sum(pos); nneg <- sum(!pos)
  if (npos == 0L || nneg == 0L || anyNA(score)) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[pos]) - npos * (npos + 1) / 2) / (npos * nneg)
}

combo_score <- function(data, features, folds) {
  prob <- rep(NA_real_, nrow(data))
  y <- as.integer(data$species == positive)
  for (idx in folds) {
    keep <- setdiff(seq_len(nrow(data)), idx)
    x_train <- as.matrix(data[keep, features, drop = FALSE])
    x_hold <- as.matrix(data[idx, features, drop = FALSE])
    fit <- xgboost::xgb.train(
      params = list(objective = "binary:logistic", max_depth = 2,
                    eta = 0.1, subsample = 0.8, colsample_bytree = 0.8,
                    nthread = 1, verbosity = 0),
      data = xgboost::xgb.DMatrix(x_train, label = y[keep]),
      nrounds = 100, verbose = 0)
    prob[idx] <- predict(fit, x_hold)
  }
  auc_rank(data$species, prob)
}

select_panel <- function(data, seed, prefix) {
  set.seed(seed)
  boruta_fit <- Boruta::Boruta(x = data[, candidate_features, drop = FALSE],
                               y = data$species, maxRuns = 500,
                               pValue = 0.05, doTrace = 0)
  boruta_fixed <- Boruta::TentativeRoughFix(boruta_fit)
  boruta_keep <- Boruta::getSelectedAttributes(boruta_fixed,
                                                withTentative = FALSE)
  if (length(boruta_keep) < 2L) stop(prefix, ": Boruta retained fewer than two predictors")

  x_lasso <- as.matrix(data[, boruta_keep, drop = FALSE])
  y_lasso <- as.integer(data$species == positive)
  set.seed(seed + 1L)
  lasso_fit <- glmnet::cv.glmnet(x_lasso, y_lasso, family = "binomial",
                                 alpha = 1, nfolds = 10,
                                 type.measure = "deviance", standardize = TRUE)
  co <- as.matrix(coef(lasso_fit, s = "lambda.min"))
  lasso_keep <- intersect(rownames(co)[co[, 1] != 0], boruta_keep)
  if (length(lasso_keep) < 2L) stop(prefix, ": LASSO retained fewer than two predictors")

  ## The AUC stage is direction-aware. The original >0.5 hard cutoff would
  ## wrongly discard inverse associations; retain the nine strongest (or all
  ## if fewer than nine) for the same combination-search stage.
  univariate <- vapply(lasso_keep, function(nm) {
    a <- auc_rank(data$species, data[[nm]])
    max(a, 1 - a)
  }, numeric(1))
  univariate <- sort(univariate, decreasing = TRUE)
  auc_keep <- names(univariate)[seq_len(min(9L, length(univariate)))]
  set.seed(seed + 2L)
  three_folds <- caret::createFolds(data$species, k = 3, returnTrain = FALSE)
  scores <- list()
  cursor <- 0L
  for (size in seq.int(2L, min(8L, length(auc_keep)))) {
    sets <- combn(auc_keep, size, simplify = FALSE)
    for (features in sets) {
      cursor <- cursor + 1L
      scores[[cursor]] <- list(features = features,
                               auc = combo_score(data, features, three_folds))
    }
  }
  values <- vapply(scores, `[[`, numeric(1), "auc")
  sizes <- vapply(scores, function(x) length(x$features), integer(1))
  order_idx <- order(-values, sizes)
  best <- scores[[order_idx[1L]]]
  if (!is.finite(best$auc)) stop(prefix, ": no finite combination AUC")

  selection <- list(panel = best$features,
                    combination_cv_auc = best$auc,
                    boruta = boruta_keep,
                    lasso = lasso_keep,
                    auc_candidates = univariate,
                    combinations_tested = length(scores),
                    boruta_fit = boruta_fixed,
                    lasso_fit = lasso_fit)
  saveRDS(selection, file.path("models", paste0(prefix, "_selection.rds")))
  write.csv(data.frame(stage = c(rep("Boruta", length(boruta_keep)),
                                   rep("LASSO", length(lasso_keep)),
                                   rep("AUC candidate", length(auc_keep)),
                                   rep("final panel", length(best$features))),
                       feature = c(boruta_keep, lasso_keep, auc_keep, best$features)),
            file.path("results", paste0(prefix, "_selection.csv")), row.names = FALSE)
  selection
}

make_recipe <- function(data) {
  recipes::recipe(species ~ ., data = data) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    themis::step_smote(species, over_ratio = 1, neighbors = 5)
}

grid_for <- function(method, p) {
  if (method == "glm") return(data.frame(parameter = "none"))
  if (method == "glmnet")
    return(expand.grid(alpha = c(0, 0.5, 1), lambda = c(0.01, 0.1)))
  if (method == "rpart") return(data.frame(cp = c(0.005, 0.02, 0.05)))
  if (method == "ranger")
    return(expand.grid(mtry = unique(pmax(1L, c(floor(sqrt(p)), floor(p/2)))),
                       splitrule = "gini", min.node.size = c(5, 15)))
  if (method == "xgbTree")
    return(expand.grid(nrounds = c(100, 250), max_depth = c(2, 3),
                       eta = c(0.05, 0.1), gamma = 0,
                       colsample_bytree = 0.8, min_child_weight = 1,
                       subsample = 0.8))
  if (method == "svmRadial")
    return(expand.grid(sigma = 1 / p, C = c(0.5, 2, 8)))
  if (method == "gbm")
    return(expand.grid(n.trees = c(100, 250),
                       interaction.depth = c(1, 3),
                       shrinkage = 0.05, n.minobsinnode = 10))
  if (method == "nnet")
    return(expand.grid(size = c(3, 5), decay = c(0.01, 0.1)))
  if (method == "knn") return(data.frame(k = c(5, 11, 21)))
  stop("Unknown model: ", method)
}

fit_method <- function(data, features, method, seed) {
  relevant <- unique(c("species", features))
  input <- data[, relevant, drop = FALSE]
  set.seed(seed)
  inner <- caret::createFolds(input$species, k = 3, returnTrain = TRUE)
  control <- caret::trainControl(method = "cv", index = inner,
                                 summaryFunction = caret::twoClassSummary,
                                 classProbs = TRUE,
                                 savePredictions = "final",
                                 allowParallel = FALSE)
  args <- list(x = make_recipe(input), data = input, method = method,
               metric = "ROC", tuneGrid = grid_for(method, length(features)),
               trControl = control)
  if (method == "glm") args$family <- binomial()
  if (method == "gbm") args$verbose <- FALSE
  if (method == "xgbTree") { args$verbose <- 0; args$nthread <- 1 }
  if (method == "nnet") { args$trace <- FALSE; args$maxit <- 200 }
  if (method == "ranger") args$num.trees <- 300
  set.seed(seed + 1L)
  do.call(caret::train, args)
}

predict_tb <- function(fit, new_data) {
  as.numeric(predict(fit, newdata = new_data, type = "prob")[[positive]])
}

## Five outer folds evaluate the *whole* feature-selection/training chain.
set.seed(2026)
outer_folds <- caret::createFolds(train_data$species, k = 5, returnTrain = FALSE)
outer_predictions <- list()
outer_auc <- list()
for (i in seq_along(outer_folds)) {
  checkpoint <- file.path("models", sprintf("outer_fold_%d.rds", i))
  if (file.exists(checkpoint)) {
    fold_result <- readRDS(checkpoint)
    message("Loaded checkpoint: outer fold ", i)
  } else {
    hold_idx <- outer_folds[[i]]
    analysis_data <- train_data[-hold_idx, , drop = FALSE]
    hold_data <- train_data[hold_idx, , drop = FALSE]
    message("Outer fold ", i, ": feature selection")
    selected <- select_panel(analysis_data, 3000L + i * 100L,
                             sprintf("outer_fold_%d", i))
    model_out <- list()
    for (j in seq_along(algorithms)) {
      method <- algorithms[j]
      message("Outer fold ", i, ": ", method,
              " on ", paste(selected$panel, collapse = ", "))
      fit <- fit_method(analysis_data, selected$panel, method,
                        4000L + i * 100L + j)
      prob <- predict_tb(fit, hold_data)
      model_out[[method]] <- list(prob = prob,
                                  auc = auc_rank(hold_data$species, prob),
                                  bestTune = fit$bestTune)
      rm(fit); gc(verbose = FALSE)
    }
    fold_result <- list(row_id = hold_data$row_id,
                        diagnosis = as.character(hold_data$species),
                        panel = selected$panel, models = model_out)
    saveRDS(fold_result, checkpoint)
  }
  for (method in algorithms) {
    item <- fold_result$models[[method]]
    outer_predictions[[length(outer_predictions) + 1L]] <- data.frame(
      fold = i, row_id = fold_result$row_id,
      diagnosis = fold_result$diagnosis, algorithm = method,
      probability_TB = item$prob)
    outer_auc[[length(outer_auc) + 1L]] <- data.frame(
      fold = i, algorithm = method, auc = item$auc,
      features = paste(fold_result$panel, collapse = ";"))
  }
}
outer_predictions <- do.call(rbind, outer_predictions)
outer_auc <- do.call(rbind, outer_auc)
write.csv(outer_predictions, "results/outer_predictions.csv", row.names = FALSE)
write.csv(outer_auc, "results/outer_auc.csv", row.names = FALSE)
ranking <- aggregate(auc ~ algorithm, outer_auc, mean)
names(ranking)[2L] <- "mean_outer_auc"
ranking <- ranking[order(-ranking$mean_outer_auc, ranking$algorithm), ]
row.names(ranking) <- NULL
write.csv(ranking, "results/model_ranking.csv", row.names = FALSE)
winner <- ranking$algorithm[1L]
message("Training-internal winner: ", winner)

## Re-run the same selector on all original training cases, then fit each
## algorithm once. Validation remains outside every fitted model; test is
## used only after the model and thresholds are locked.
full_selection <- select_panel(train_data, seed = 9991L, prefix = "full_training")
panel <- full_selection$panel
final_models <- list()
validation_predictions <- list()
test_predictions <- list()
for (j in seq_along(algorithms)) {
  method <- algorithms[j]
  message("Final training model: ", method)
  fit <- fit_method(train_data, panel, method, 10000L + j)
  final_models[[method]] <- fit
  saveRDS(fit, file.path("models", paste0("final_", method, ".rds")))
  validation_predictions[[method]] <- data.frame(
    row_id = prep$validation$row_id,
    diagnosis = as.character(prep$validation$species),
    algorithm = method,
    probability_TB = predict_tb(fit, prep$validation))
  test_predictions[[method]] <- data.frame(
    row_id = prep$test$row_id,
    diagnosis = as.character(prep$test$species),
    algorithm = method,
    probability_TB = predict_tb(fit, prep$test))
}
validation_predictions <- do.call(rbind, validation_predictions)
test_predictions <- do.call(rbind, test_predictions)
write.csv(validation_predictions, "results/validation_predictions.csv", row.names = FALSE)
write.csv(test_predictions, "results/test_predictions.csv", row.names = FALSE)

## Required demographic comparators, using the winning algorithm on the
## same frozen allocation and the same resampling/preprocessing policy.
comparators <- list(age_sex = c("Age", "Gender"),
                    panel_age_sex = c(panel, "Age", "Gender"))
for (j in seq_along(comparators)) {
  name <- names(comparators)[j]
  message("Demographic comparator: ", name)
  fit <- fit_method(train_data, comparators[[j]], winner, 12000L + j)
  saveRDS(fit, file.path("models", paste0(name, "_", winner, ".rds")))
  validation_predictions <- rbind(validation_predictions, data.frame(
    row_id = prep$validation$row_id,
    diagnosis = as.character(prep$validation$species), algorithm = name,
    probability_TB = predict_tb(fit, prep$validation)))
  test_predictions <- rbind(test_predictions, data.frame(
    row_id = prep$test$row_id,
    diagnosis = as.character(prep$test$species), algorithm = name,
    probability_TB = predict_tb(fit, prep$test)))
}

## The historical five-feature XGBoost is evaluated even if the reselected
## panel or winning algorithm changes. It is a prespecified comparator, not
## substituted for the actually selected main analysis.
message("Historical five-feature XGBoost comparator")
legacy_fit <- fit_method(train_data, old_panel, "xgbTree", 13001L)
saveRDS(legacy_fit, "models/historical_five_xgbTree.rds")
validation_predictions <- rbind(validation_predictions, data.frame(
  row_id = prep$validation$row_id,
  diagnosis = as.character(prep$validation$species),
  algorithm = "historical_five_xgbTree",
  probability_TB = predict_tb(legacy_fit, prep$validation)))
test_predictions <- rbind(test_predictions, data.frame(
  row_id = prep$test$row_id,
  diagnosis = as.character(prep$test$species),
  algorithm = "historical_five_xgbTree",
  probability_TB = predict_tb(legacy_fit, prep$test)))
write.csv(validation_predictions, "results/validation_predictions.csv", row.names = FALSE)
write.csv(test_predictions, "results/test_predictions.csv", row.names = FALSE)

saveRDS(list(winner = winner, panel = panel, legacy_panel = old_panel,
             candidate_features = candidate_features,
             outer_auc = outer_auc, model_ranking = ranking,
             training_n = nrow(train_data),
             validation_n = nrow(prep$validation),
             test_n = nrow(prep$test)), "results/analysis_lock.rds")
writeLines(c(sprintf("Training-internal winning algorithm: %s", winner),
             sprintf("Final selected panel: %s", paste(panel, collapse = ", ")),
             sprintf("Historical panel: %s", paste(old_panel, collapse = ", ")),
             "Five outer folds; Boruta maxRuns=500; LASSO 10 folds;",
             "direction-aware univariate AUC top nine; XGBoost-scored 2-8-feature search with 3 folds;",
             "nine algorithms with three-fold internal hyperparameter tuning;",
             "SMOTE inside each training fold only; TB is the positive class."),
           "results/model_lock.txt", useBytes = TRUE)
print(ranking)
cat("Selected panel:", paste(panel, collapse = ", "), "\n")
