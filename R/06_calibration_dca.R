# ============================================================
# 06_calibration_dca.R
# Calibration curve and Decision Curve Analysis
# Corresponds to: Figure 4A (calibration), 4B (DCA)
# ============================================================

library(rms)
library(dcurves)
library(ggplot2)
library(dplyr)

# ---- Calibration Curve (Figure 4A) ----
actual_y_val <- ifelse(y_val_raw == "TB", 1, 0)

pdf("output/Figure4A_Calibration.pdf", width = 7, height = 7)
val.prob(p = val_pred, y = actual_y_val, pl = TRUE, smooth = TRUE,
         logistic.cal = FALSE,
         statloc = c(0.0, 1), cex = 0.8,
         xlab = "Predicted Probability of TB",
         ylab = "Actual Observed Fraction of TB")
title(main = "Calibration Curve (Validation Set)", font.main = 2, cex.main = 1.2)
dev.off()

# ---- Decision Curve Analysis (Figure 4B) ----
dca_data <- data.frame(
  tb_status = actual_y_val,
  xgb_prob = val_pred
)

dca_result <- dca(tb_status ~ xgb_prob, data = dca_data,
                  thresholds = seq(0, 0.8, by = 0.01))

p_dca <- plot(dca_result) +
  theme_classic() +
  scale_color_manual(
    values = c("xgb_prob" = "#E64B35", "Treat All" = "#95A5A6",
               "Treat None" = "black"),
    labels = c("xgb_prob" = "XGBoost", "Treat All" = "Treat All as TB",
               "Treat None" = "Treat None")
  ) +
  labs(title = "Decision Curve Analysis (Validation Set)",
       x = "Threshold Probability", y = "Net Benefit") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 15),
        axis.text = element_text(size = 14, color = "black"),
        legend.position = c(0.75, 0.75))

ggsave("output/Figure4B_DCA.pdf", p_dca, width = 8, height = 6)
