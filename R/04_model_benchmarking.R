# ============================================================
# 04_model_benchmarking.R
# 9-model comparison + exhaustive combination search
# Corresponds to: Figure 3C (CV boxplot)
# ============================================================

library(caret)
library(doParallel)
library(pROC)
library(ggplot2)
library(tidyr)
library(dplyr)

selected_cols <- c("species", valid_vars)
train_elite <- train_s[, selected_cols]
val_elite   <- val_s[, selected_cols]

train_elite$species <- factor(train_elite$species,
                               levels = levels(TB_Pneumonia$species))
val_elite$species   <- factor(val_elite$species,
                               levels = levels(TB_Pneumonia$species))

x_train <- train_elite[, setdiff(names(train_elite), "species"), drop = FALSE]
y_train <- train_elite$species

# ---- 5-fold CV control ----
ctrl <- trainControl(method = "cv", number = 5,
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE, verboseIter = FALSE,
                     allowParallel = TRUE)

# ---- 9 algorithms ----
model_methods <- c(glm = "glm", glmnet = "glmnet", rpart = "rpart",
                   ranger = "ranger", xgbTree = "xgbTree",
                   svmRadial = "svmRadial", gbm = "gbm",
                   nnet = "nnet", knn = "knn")
tune_len <- list(glm = 1, glmnet = 5, rpart = 5, ranger = 5,
                 xgbTree = 5, svmRadial = 5, gbm = 5, nnet = 5, knn = 5)

cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

fits <- list()
for (nm in names(model_methods)) {
  method <- model_methods[[nm]]
  extra_args <- list()
  if (method == "glm") extra_args$family <- binomial()
  if (method == "xgbTree") { extra_args$verbose <- 0; extra_args$verbosity <- 0 }
  if (method == "gbm") extra_args$verbose <- FALSE
  if (method == "nnet") extra_args$trace <- FALSE

  set.seed(2026)
  args <- c(list(x = x_train, y = y_train, method = method,
                 metric = "ROC", trControl = ctrl,
                 tuneLength = tune_len[[nm]]), extra_args)
  fit <- try(suppressWarnings(do.call(caret::train, args)), silent = TRUE)
  if (!inherits(fit, "try-error")) fits[[nm]] <- fit
}
stopCluster(cl)
registerDoSEQ()

# ---- CV Performance Boxplot (Figure 3C) ----
resamps <- resamples(fits)
plot_data <- resamps$values
plot_long <- pivot_longer(plot_data, cols = -Resample,
                          names_to = "Model_Metric", values_to = "Value")
plot_roc <- plot_long[grepl("~ROC$", plot_long$Model_Metric), ]
plot_roc$Model <- gsub("~ROC", "", plot_roc$Model_Metric)

model_order <- plot_roc %>%
  group_by(Model) %>% summarise(med = median(Value, na.rm = TRUE)) %>%
  arrange(med) %>% pull(Model)
plot_roc$Model <- factor(plot_roc$Model, levels = model_order)

top3 <- tail(model_order, 3)
plot_roc$FillColor <- ifelse(plot_roc$Model %in% top3, "#E64B35", "#7CB5EC")

p_cv <- ggplot(plot_roc, aes(x = Value, y = Model)) +
  stat_boxplot(geom = "errorbar", width = 0.3, color = "#2B5B84") +
  geom_boxplot(aes(fill = FillColor), color = "#2B5B84",
               outlier.shape = 1, outlier.color = "#2B5B84") +
  stat_summary(fun = median, geom = "point", shape = 16,
               size = 4, color = "black") +
  scale_fill_identity() +
  labs(title = "Cross-Validation ROC Across 9 ML Algorithms",
       x = "ROC", y = "") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.text = element_text(size = 15, color = "black"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

ggsave("output/Figure3C_CV_Boxplot.pdf", p_cv, width = 10, height = 6)

# ---- Exhaustive combination search (2-8 features) ----
cat("Running exhaustive combination search...\n")
max_combo <- min(8, length(valid_vars))
all_combinations <- unlist(lapply(2:max_combo,
  function(i) combn(valid_vars, m = i, simplify = FALSE)), recursive = FALSE)

fast_ctrl <- trainControl(method = "cv", number = 3,
                          summaryFunction = twoClassSummary,
                          classProbs = TRUE, allowParallel = TRUE)
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

combo_results <- data.frame(Combo_ID = integer(), Num_Features = integer(),
                            Features = character(), Val_AUC = numeric(),
                            stringsAsFactors = FALSE)

for (idx in seq_along(all_combinations)) {
  current_features <- all_combinations[[idx]]
  x_train_sub <- train_s[, current_features, drop = FALSE]
  x_val_sub   <- val_s[, current_features, drop = FALSE]
  set.seed(2026)
  fit_xgb <- try(caret::train(x = x_train_sub, y = train_s$species,
                method = "xgbTree", metric = "ROC",
                trControl = fast_ctrl, tuneLength = 3,
                verbose = 0, verbosity = 0), silent = TRUE)
  if (!inherits(fit_xgb, "try-error")) {
    pr <- predict(fit_xgb, newdata = x_val_sub, type = "prob")
    roc_val <- pROC::roc(val_s$species, pr[["TB"]],
                          levels = levels(val_s$species),
                          direction = "<", quiet = TRUE)
    combo_results <- rbind(combo_results, data.frame(
      Combo_ID = idx, Num_Features = length(current_features),
      Features = paste(current_features, collapse = " + "),
      Val_AUC = as.numeric(pROC::auc(roc_val))))
  }
}
stopCluster(cl); registerDoSEQ()

combo_results <- combo_results[order(-combo_results$Val_AUC), ]
write.csv(combo_results, "output/All_Combinations_Results.csv", row.names = FALSE)
cat("Best combo:", combo_results$Features[1],
    "AUC:", combo_results$Val_AUC[1], "\n")
