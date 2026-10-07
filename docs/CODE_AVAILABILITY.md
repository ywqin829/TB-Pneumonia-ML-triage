# Code Availability Statement

The code used for model training, validation, and visualization in this study is publicly available in the following repository:

> **GitHub**: [https://github.com/ywqin829/TB-Pneumonia-ML-triage](https://github.com/ywqin829/TB-Pneumonia-ML-triage)

The repository is archived on Zenodo (concept DOI, always resolves to the latest version):

> **Zenodo**: [https://doi.org/10.5281/zenodo.21151460](https://doi.org/10.5281/zenodo.21151460)

## Contents

| File | Description |
|---|---|
| `R/00_prepare.R` | Source audit, derivable index reconstruction, and frozen diagnosis-stratified 60/20/20 split (seed 123) |
| `R/01_model_pipeline.R` | Boruta and LASSO feature engineering, direction-aware AUC screening, exhaustive feature-combination search, and nine-algorithm benchmarking with all learned steps inside training folds |
| `R/02_evaluate.R` | Locked validation-derived thresholds and one-pass test-set evaluation |
| `R/03_supplement.R` | Supplemental descriptive, stability, demographic, and reproducibility outputs |
| `R/04_figures.R` | Scientific figures rebuilt from the locked outputs (exported per panel) |
| `R/05_integrity_check.R` | Independent, read-only audit of saved revision outputs |
| `R/05_smoke_check.R` | Delivery smoke check of saved predictions, unit conversion, and the Shiny server |
| `R/06_balanced_evaluation.R` | Author-selected balanced operating strategy (disclosed post-result adjustment) |
| `R/app.R` | Interactive Shiny research dashboard |
| `R/app_predict.R` | Prediction utilities for the dashboard |
| `data/` | Directory for source data (not distributed due to privacy restrictions) |

## Dependencies

All analyses were conducted in **R 4.5.1** with the following key packages:

- `rsample`, `caret`, `recipes`, `themis` — data partitioning, model training, and cross-validation
- `xgboost`, `ranger`, `glmnet`, `gbm`, `nnet`, `rpart`, `kernlab` — candidate algorithms
- `Boruta` — feature selection
- `pROC` — ROC curve analysis and threshold selection
- `ggplot2` — figure generation
- `shiny` — interactive dashboard

## Usage

1. Place the source dataset as `TB_Pneumonia.csv` in the working directory (contact the corresponding author for access).
2. Run the scripts in order from `00_prepare.R` through `06_balanced_evaluation.R`.
3. The dashboard (`app.R`) requires the saved `results/` and `models/` outputs produced by the pipeline.

## Interactive Web App

A live version of the research decision-support prototype is available at:

> [https://nana2379723224.shinyapps.io/tb-pneumonia-research/](https://nana2379723224.shinyapps.io/tb-pneumonia-research/)

## License

This code is made available for research and educational purposes. For clinical use inquiries, please contact the corresponding author.
