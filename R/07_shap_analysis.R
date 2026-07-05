# ============================================================
# 07_shap_analysis.R
# SHAP beeswarm and individual force plots
# Corresponds to: Figure 5A (beeswarm), 5B (force plots)
# ============================================================

library(shapviz)
library(ggplot2)
library(patchwork)

# Extract raw xgboost model
xgb_raw <- xgb_final$finalModel

# ---- Global SHAP Beeswarm (Figure 5A) ----
X_test_mat <- as.matrix(x_test_raw)
sv_global <- shapviz(xgb_raw, X_pred = X_test_mat, X = x_test_raw)

# Flip sign: positive SHAP = TB
sv_global$S <- -sv_global$S
sv_global$baseline <- -sv_global$baseline

p_beeswarm <- sv_importance(sv_global, kind = "beeswarm") +
  scale_color_viridis_c(option = "magma", begin = 0.1, end = 0.9) +
  theme_bw() +
  labs(title = "SHAP Summary (Impact on TB Probability)",
       x = "SHAP Value") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text.y = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(color = "black", size = 11),
        legend.position = "bottom")

ggsave("output/Figure5A_SHAP_Beeswarm.pdf", p_beeswarm, width = 8, height = 5)

# ---- Individual SHAP Force Plots (Figure 5B) ----
X_val_mat <- as.matrix(x_val_raw)
id_high_tb   <- which.max(val_pred)
id_high_pneu <- which.min(val_pred)

sv_high <- shapviz(xgb_raw,
                   X_pred = X_val_mat[id_high_tb, , drop = FALSE],
                   X = x_val_raw[id_high_tb, , drop = FALSE])
sv_low  <- shapviz(xgb_raw,
                   X_pred = X_val_mat[id_high_pneu, , drop = FALSE],
                   X = x_val_raw[id_high_pneu, , drop = FALSE])

sv_high$S <- -sv_high$S; sv_high$baseline <- -sv_high$baseline
sv_low$S  <- -sv_low$S;  sv_low$baseline  <- -sv_low$baseline

prob_A <- round(val_pred[id_high_tb] * 100, 1)
prob_B <- round(val_pred[id_high_pneu] * 100, 1)

p_force_tb <- sv_force(sv_high) + theme_minimal() +
  labs(title = paste0("Patient A: High TB Risk (", prob_A, "%)")) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5,
                                  color = "#E64B35"))

p_force_pneu <- sv_force(sv_low) + theme_minimal() +
  labs(title = paste0("Patient B: Low TB Risk (", prob_B, "%)"),
       x = "SHAP Value") +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5,
                                  color = "#2980B9"))

ggsave("output/Figure5B_Force_Plots.pdf",
       p_force_tb / p_force_pneu, width = 10, height = 6)
