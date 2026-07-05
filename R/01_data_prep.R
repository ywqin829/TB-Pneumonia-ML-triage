# ============================================================
# 01_data_prep.R
# Data Loading, Cleaning, and Table 1 (Baseline Characteristics)
# Corresponds to: Supplementary Table S1, Figure 2A-E data prep
# ============================================================

# ---- Package loading ----
if (!requireNamespace("xfun", quietly = TRUE)) install.packages("xfun")
pkgs <- c("tidyverse", "tableone", "janitor", "skimr")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Data loading ----
# Update the path to your local CSV file
base_path <- "data/TB_Pneumonia_with_new_features.csv"
TB_Pneumonia <- read_csv(base_path, show_col_types = FALSE) %>%
  filter(complete.cases(.))

# ---- Cleaning ----
TB_Pneumonia$species <- as.factor(TB_Pneumonia$species)
TB_Pneumonia <- remove_empty(TB_Pneumonia, which = c("rows", "cols"))
if (sum(duplicated(TB_Pneumonia)) > 0) {
  TB_Pneumonia <- TB_Pneumonia[!duplicated(TB_Pneumonia), ]
}

cat("Final cohort size:", nrow(TB_Pneumonia), "\n")
print(table(TB_Pneumonia$species))

# ---- Table 1: Baseline Characteristics ----
all_vars <- setdiff(names(TB_Pneumonia), "species")

# Define non-normal variables for median [IQR] reporting
non_normal_vars <- c("NLR", "PLR", "NPR", "dNLR", "MLR", "NMLR",
                     "SIRI", "SII", "LMR", "MNR", "PNR", "IPIV", "PIV", "AIS")

table_one <- CreateTableOne(
  vars = all_vars,
  strata = "species",
  data = TB_Pneumonia
)

table1 <- print(table_one,
                nonnormal = non_normal_vars,
                showAllLevels = TRUE,
                noSpaces = TRUE,
                printToggle = TRUE,
                quote = FALSE,
                exact = FALSE,
                smd = TRUE,
                explain = TRUE)

write.csv(as.data.frame(table1), "output/Table_1_Baseline_Characteristics.csv")
