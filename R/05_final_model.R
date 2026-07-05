# ============================================================
# 05_final_model.R
# Final XGBoost model training, ROC curves, confusion matrix
# Corresponds to: Figure 3D (ROC curves), 3E (confusion matrix)
# ============================================================

library(caret)
library(xgboost)
library(pROC)
library(ggplot2)
library(dplyr)

# ---- Load enhanced training data ----
TB_NEW <- read.csv("data/TB_NEW1.csv", header = TRUE) %>% na.omit()
final_vars <- c("RBC", "PDW", "MONO.", "MPV", "LMR")

train_old_clean <- train_raw %>%
  select(all_of(final_vars), species)
train_enhanced <- rbind(train_old_clean, TB_NEW)
cat("Enhanced training set:", nrow(train_enhanced), "samples\n")

# ---- XGBoost hyperparameter grid ----
xgb_grid <- expand.grid(
  nrounds = c(100, 150, 200),
  max_depth = c(3, 5),
  eta = c(0.01, 0.05),
  gamma = 0,
  colsample_bytree = c(0.8, 1.0),
  min_child_weight = c(1, 3),
  subsample = c(0.8, 1.0)
)

ctrl <- trainControl(method = "cv", number = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(123)
xgb_final <- train(
  species ~ ., data = train_enhanced,
  method = "xgbTree", trControl = ctrl,
  tuneGrid = xgb_grid, metric = "ROC"
)

print(xgb_final$bestTune)
cat("Best CV AUC:", max(xgb_final$results$ROC), "\n")

# ---- ROC Curves (Figure 3D) ----
x_val_raw <- val_raw[, final_vars, drop = FALSE]
y_val_raw <- factor(val_raw$species, levels = levels(train_raw$species))
x_test_raw <- test_raw[, final_vars, drop = FALSE]
y_test_raw <- factor(test_raw$species, levels = levels(train_raw$species))

train_pred <- predict(xgb_final, newdata = train_enhanced[, final_vars],
                      type = "prob")[["TB"]]
val_pred   <- predict(xgb_final, newdata = x_val_raw, type = "prob")[["TB"]]
test_pred  <- predict(xgb_final, newdata = x_test_raw, type = "prob")[["TB"]]

roc_train <- pROC::roc(train_enhanced$species, train_pred,
                        levels = levels(y_val_raw), direction = "<", quiet = TRUE)
roc_val   <- pROC::roc(y_val_raw, val_pred,
                        levels = levels(y_val_raw), direction = "<", quiet = TRUE)
roc_test  <- pROC::roc(y_test_raw, test_pred,
                        levels = levels(y_val_raw), direction = "<", quiet = TRUE)

p_roc <- ggroc(list(Train = roc_train, Validation = roc_val, Test = roc_test),
               linewidth = 1, legacy.axes = TRUE) +
  theme_classic() +
  scale_color_manual(values = c(Train = "#95A5A6", Validation = "#2980B9",
                                 Test = "#E64B35")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "ROC Curves", x = "1 - Specificity", y = "Sensitivity") +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, color = "black"),
        legend.position = c(0.50, 0.25))

ggsave("output/Figure3D_ROC.pdf", p_roc, width = 8, height = 7)

# ---- Confusion Matrix (Figure 3E) ----
optimal_coords <- pROC::coords(roc_val, "best",
                                best.method = "youden",
                                ret = c("threshold", "specificity",
                                        "sensitivity", "accuracy"))

pred_class <- ifelse(test_pred >= optimal_coords$threshold, "TB", "Pneumonia")
pred_class <- factor(pred_class, levels = levels(y_test_raw))
cm <- confusionMatrix(pred_class, y_test_raw, positive = "TB")
cm_data <- as.data.frame(cm$table)

p_cm <- ggplot(cm_data, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = Freq), fontface = "bold", size = 10) +
  scale_fill_gradient(low = "#F0F8FF", high = "#4DBBD5") +
  theme_minimal() +
  labs(title = "Confusion Matrix (Test Set)",
       subtitle = sprintf("Cut-off: %.3f | Sens: %.1f%% | Spec: %.1f%%",
                          optimal_coords$threshold,
                          optimal_coords$sensitivity * 100,
                          optimal_coords$specificity * 100)) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 16, color = "black"),
        panel.grid = element_blank())

ggsave("output/Figure3E_Confusion_Matrix.pdf", p_cm, width = 6, height = 5)
