# ============================================================
# 02_exploratory_figures.R
# Baseline demographic plots and immune-inflammatory profiles
# Corresponds to: Figure 2A (Age pyramid), 2B (Gender bar),
#                 2C (Radar chart), 2D (Correlation heatmap),
#                 2E (Boruta plot)
# ============================================================

library(ggplot2)
library(dplyr)
library(corrplot)
library(fmsb)

# ---- Figure 2A: Age Distribution Pyramid ----
TB_Pneumonia$AgeGroup <- cut(TB_Pneumonia$Age, breaks = seq(0, 100, by = 5))

pyramid_data <- TB_Pneumonia %>%
  group_by(species, AgeGroup) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(
    percentage = n / sum(n) * 100,
    percentage_label = paste0(round(percentage, 1), "%")
  ) %>%
  ungroup() %>%
  mutate(percentage_plot = ifelse(species == "TB", -percentage, percentage))

p_pyramid <- ggplot(pyramid_data,
                    aes(x = AgeGroup, y = percentage_plot, fill = species)) +
  geom_bar(stat = "identity", width = 0.9, color = "white") +
  geom_line(aes(y = percentage_plot, group = species),
            color = "black", linetype = "dashed", linewidth = 0.3) +
  geom_text(aes(label = ifelse(percentage > 0.5, percentage_label, ""),
                hjust = ifelse(species == "TB", 1.1, -0.1)),
            size = 3.0, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("TB" = "#E64B35", "Pneumonia" = "#4DBBD5")) +
  scale_y_continuous(limits = c(-20, 20),
                     breaks = seq(-20, 20, 5),
                     labels = function(x) paste0(abs(x), "%")) +
  labs(title = "Age Distribution Pyramid", x = "Age Group",
       y = "Percentage within Each Disease") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.text = element_text(size = 11), legend.position = "bottom")

ggsave("output/Figure2A_Age_Pyramid.pdf", p_pyramid, width = 8, height = 6)

# ---- Figure 2B: Gender Distribution ----
gender_data <- TB_Pneumonia %>%
  mutate(Gender = factor(Gender, levels = c(1, 2), labels = c("Male", "Female"))) %>%
  group_by(species, Gender) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(percentage = n / sum(n) * 100,
         percentage_label = paste0(round(percentage, 1), "%")) %>%
  ungroup()

p_gender <- ggplot(gender_data,
                   aes(x = percentage, y = species, fill = Gender)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.5) +
  geom_text(aes(label = percentage_label),
            position = position_dodge(width = 0.7),
            hjust = -0.2, size = 4, fontface = "bold") +
  scale_x_continuous(limits = c(0, max(gender_data$percentage) + 15),
                     labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = c("Male" = "#3C5488", "Female" = "#F39B7F")) +
  labs(title = "Gender Distribution") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text = element_text(size = 14, color = "black", face = "bold"),
        legend.position = "top", legend.title = element_blank())

ggsave("output/Figure2B_Gender.pdf", p_gender, width = 8, height = 4)

# ---- Figure 2C: Radar Chart of Inflammatory Indices ----
radar_vars <- c("NLR", "PLR", "LMR", "SII", "SIRI", "AIS", "dNLR", "MLR")
radar_df <- TB_Pneumonia %>%
  group_by(species) %>%
  summarise(across(all_of(radar_vars), median), .groups = "drop") %>%
  column_to_rownames("species")

radar_scaled <- as.data.frame(
  apply(radar_df, 2, function(x) (x - min(x)) / (max(x) - min(x) + 1e-8))
)
radar_plot <- rbind(rep(1, ncol(radar_scaled)),
                    rep(0, ncol(radar_scaled)),
                    radar_scaled)
rownames(radar_plot) <- c("max", "min", "TB", "Pneumonia")

pdf("output/Figure2C_Radar.pdf", width = 7, height = 7)
par(mar = c(1, 1, 2, 1))
radarchart(radar_plot, axistype = 1,
           pcol = c("#C0392B", "#2980B9"),
           pfcol = scales::alpha(c("#C0392B", "#2980B9"), 0.3),
           plwd = 2.5, plty = 1, cglcol = "grey80", cglty = 2,
           axislabcol = "grey30", vlcex = 1.2,
           title = "Inflammatory Index Radar — Median Profiles")
legend("topright", legend = c("TB", "Pneumonia"),
       col = c("#E64B35", "#4DBBD5"), lty = 1, lwd = 2.5, bty = "n", cex = 1.1)
dev.off()

# ---- Figure 2D: Spearman Correlation Heatmap ----
vars_for_cor <- setdiff(names(TB_Pneumonia), c("species", "Gender"))
cor_tb   <- cor(TB_Pneumonia[TB_Pneumonia$species == "TB", vars_for_cor],
                method = "spearman", use = "complete.obs")
cor_pneu <- cor(TB_Pneumonia[TB_Pneumonia$species == "Pneumonia", vars_for_cor],
                method = "spearman", use = "complete.obs")
cor_pooled <- cor(TB_Pneumonia[, vars_for_cor],
                  method = "spearman", use = "complete.obs")
hc_order <- hclust(as.dist(1 - cor_pooled))$order

pdf("output/Figure2D_Heatmap.pdf", width = 14, height = 7)
par(mfrow = c(1, 2), mar = c(1, 1, 2, 1), oma = c(0, 0, 3, 0))
corrplot(cor_tb[hc_order, hc_order],
         method = "color", type = "upper", order = "original",
         col = COL2("RdBu", 200), col.lim = c(-1, 1),
         tl.cex = 0.65, tl.col = "black", title = "TB Group")
corrplot(cor_pneu[hc_order, hc_order],
         method = "color", type = "upper", order = "original",
         col = COL2("RdBu", 200), col.lim = c(-1, 1),
         tl.cex = 0.65, tl.col = "black", title = "Pneumonia Group")
mtext("Spearman Correlation Heatmap", side = 3, line = 1,
      outer = TRUE, cex = 1.2, font = 2)
dev.off()
par(mfrow = c(1, 1))
