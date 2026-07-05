# Code Availability Statement

The code used for model training, validation, and visualization in this study is publicly available in the following repository:

> **GitHub**: [https://github.com/ywqin829/TB-Pneumonia-ML-Diagnosis-Triage](https://github.com/ywqin829/TB-Pneumonia-ML-Diagnosis-Triage)

## Contents

| Directory / File | Description |
|---|---|
| `R/01_data_prep.R` | Data loading, cleaning, calculation of derived indices, and Table 1 generation |
| `R/02_exploratory_figures.R` | Exploratory visualizations: age pyramid, gender distribution, radar plot, correlation heatmap (Figure 2) |
| `R/03_feature_selection.R` | Boruta, LASSO regression, and univariate AUC screening (Figures 2E, 3A, 3B) |
| `R/04_model_benchmarking.R` | 9-model comparison with 5-fold CV and exhaustive feature combination search (Figure 3C) |
| `R/05_final_model.R` | Final XGBoost training with 10-fold CV, ROC curves, and confusion matrix (Figures 3D, 3E) |
| `R/06_calibration_dca.R` | Calibration curve and decision curve analysis (Figures 4A, 4B) |
| `R/07_shap_analysis.R` | SHAP beeswarm summary plot and individual force plots (Figures 5A, 5B) |
| `R/08_shiny_app.R` | Interactive Shiny dashboard for clinical decision support |
| `data/` | Directory for raw data (not publicly distributed due to privacy restrictions) |
| `output/` | Directory for generated figures and tables |

## Dependencies

All analyses were conducted in **R 4.x** with the following key packages:

- `caret`, `xgboost` — model training and cross-validation
- `pROC` — ROC curve analysis
- `glmnet` — LASSO regression
- `Boruta` — feature selection
- `shapviz` — SHAP explainability
- `dcurves` — decision curve analysis
- `shiny`, `bslib`, `echarts4r` — interactive dashboard

## Usage

1. Place the raw dataset in the `data/` directory.
2. Source scripts sequentially from `01_data_prep.R` through `07_shap_analysis.R`.
3. The Shiny app (`08_shiny_app.R`) requires the pre-trained model file `tb_pneumonia_model_V2.rds`.

## Interactive Web App

A live version of the diagnostic dashboard is available at:

> [https://0pwfel-0-0.shinyapps.io/TB_AI_App/](https://0pwfel-0-0.shinyapps.io/TB_AI_App/)

## License

This code is made available for research and educational purposes. For clinical use inquiries, please contact the corresponding author.
