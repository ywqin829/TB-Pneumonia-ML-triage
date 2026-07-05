# ============================================================
# 03_feature_selection.R
# Boruta + LASSO + Univariate AUC screening
# Corresponds to: Figure 2E (Boruta), Figure 3A (LASSO path),
#                 Figure 3B (univariate AUC lollipop)
# ============================================================

library(dplyr)
library(Boruta)
library(glmnet)
library(pROC)
library(ggplot2)
library(caret)

# Load preprocessed data (run 01_data_prep.R first, then this)
# Assumes: train_s, val_s, test_s are available from 01_data_prep.R
if (!exists("train_s")) stop("Run 01_data_prep.R and data splitting first.")

# ---- Step 1: Boruta Feature Selection ----
set.seed(2026)
boruta_res <- Boruta(species ~ ., data = train_s,
                     pValue = 0.05, mcAdj = TRUE, maxRuns = 500)

# Boruta importance plot (Figure 2E)
pdf("output/Figure2E_Boruta.pdf", width = 10, height = 7)
par(mar = c(6, 4, 1.5, 2) + 0.1)
plot(boruta_res, sort = TRUE, xlab = "",
     main = "Boruta Feature Importance",
     las = 2, cex.axis = 1, cex.names = 0.72)
dev.off()

boruta_vars <- getSelectedAttributes(boruta_res, withTentative = FALSE)
cat("Boruta confirmed:", length(boruta_vars), "features\n")

# ---- Step 2: LASSO Regression (10-fold CV) ----
x_train <- as.matrix(train_s[, setdiff(names(train_s), "species")])
y_train <- train_s$species

set.seed(2026)
cv_lasso <- cv.glmnet(x_train, y_train, family = "binomial",
                       alpha = 1, nfolds = 10)

# LASSO path plot (Figure 3A)
pdf("output/Figure3A_LASSO.pdf", width = 8, height = 6)
plot(cv_lasso)
title("LASSO 10-Fold Cross-Validation", line = 2.5)
dev.off()

lasso_coef <- coef(cv_lasso, s = "lambda.min")
lasso_vars <- setdiff(rownames(lasso_coef)[which(lasso_coef != 0)], "(Intercept)")
cat("LASSO selected:", length(lasso_vars), "features\n")

# ---- Step 3: Univariate AUC Filtering ----
single_aucs <- list()
for (var in lasso_vars) {
  roc_obj <- pROC::roc(response = train_s$species,
                       predictor = train_s[[var]],
                       levels = levels(train_s$species),
                       direction = "<", quiet = TRUE)
  single_aucs[[var]] <- as.numeric(pROC::auc(roc_obj))
}

auc_df <- data.frame(Feature = names(single_aucs), AUC = unlist(single_aucs))
auc_df <- auc_df[order(-auc_df$AUC), ]
valid_vars <- auc_df$Feature[auc_df$AUC > 0.5]
cat("After AUC>0.5 filter:", length(valid_vars), "features\n")

# Lollipop chart (Figure 3B)
plot_df <- auc_df[auc_df$Feature %in% valid_vars, ]
plot_df$Feature <- reorder(plot_df$Feature, plot_df$AUC)

p_lollipop <- ggplot(plot_df, aes(x = AUC, y = Feature)) +
  geom_segment(aes(x = 0.5, xend = AUC, y = Feature, yend = Feature),
               color = "grey60", linewidth = 1.2) +
  geom_point(color = "#E64B35", size = 5) +
  geom_text(aes(label = sprintf("%.3f", AUC)),
            hjust = -0.4, size = 4.5, fontface = "bold") +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "grey40", linewidth = 1) +
  scale_x_continuous(limits = c(0.48, max(plot_df$AUC) + 0.05)) +
  labs(title = "Single-Variable AUC", x = "AUC", y = "") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0),
        axis.text.y = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"))

ggsave("output/Figure3B_Univariate_AUC.pdf", p_lollipop, width = 8, height = 5)
