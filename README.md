# TB-Pneumonia-ML-triage

**Machine learning-based triage of tuberculosis versus pneumonia using routine complete blood count parameters.**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21151460.svg)](https://doi.org/10.5281/zenodo.21151460)

## Overview

This repository contains the analysis code for an explainable machine learning model that distinguishes pulmonary tuberculosis (TB) from pneumonia using three routine complete blood count (CBC) parameters: MCV, PDW, and MONO%. The analysis follows the route **Boruta → LASSO → univariate AUC screening → feature-combination screening → benchmarking of nine machine learning algorithms**.

The raw source data are read from `TB_Pneumonia.csv` in the working directory (see [Requirements](#requirements)). De-identified patient data are not distributed; please contact the corresponding author for access.

## Key Results

| Item | Value |
|---|---|
| Cohort | 2,382 records (711 TB; 1,671 pneumonia) |
| Data split (6:2:2, stratified) | Training 1,428 / Validation 476 / Test 478 |
| Selected panel | MCV, PDW, MONO% |
| Model | XGBoost |
| Test AUC | 0.833 (95% CI, 0.793–0.873) |
| Primary cutoff (validation Youden) | 0.422576 |
| Sensitivity | 84.6% (95% CI, 77.8%–89.6%) |
| Specificity | 71.0% (95% CI, 66.0%–75.6%) |
| Balanced accuracy | 77.8% (95% CI, 73.9%–81.5%) |
| Brier score | 0.166 |

A supplementary dual-threshold triage strategy is retained: sensitivity ≥ 90% cutoff at 0.290 and rule-in specificity ≥ 90% cutoff at 0.634, defining lower (< 0.290), intermediate (0.290–0.634), and higher (≥ 0.634) score zones.

## Analysis Scripts

| Script | Function |
|---|---|
| `00_prepare.R` | Read the raw data; check for missing values and completely duplicate rows; assemble the candidate variables; compute LMR and PIV; and split the cohort by diagnosis into stratified random 6:2:2 training, validation, and test sets. |
| `01_model_pipeline.R` | Perform Boruta, LASSO, univariate AUC ranking, and feature-combination screening; repeat the selection independently within each of the five outer training partitions; compare nine algorithms; and save the final model together with per-dataset predictions. |
| `02_evaluate.R` | Determine the classification threshold on the validation set; compute test-set AUC, the confusion matrix, sensitivity, specificity, and confidence intervals; and carry out paired AUC comparison, calibration, decision-curve, subgroup, and score-zone analyses. |
| `03_supplement.R` | Automatically source `02_evaluate.R`; supplement baseline statistics for diagnosis groups and data splits, feature-selection stability, demographic comparisons, final hyperparameters, and bootstrap intervals for the calibration parameters. |
| `04_figures.R` | Compute SHAP contributions from the saved data, models, and results; draw the study flow, baseline, feature-selection, ROC, calibration, decision-curve, and SHAP figures; and export each panel separately as PDF and PNG. |
| `05_integrity_check.R` | Independently verify the consistency of the saved data splits, predictions, performance metrics, and statistical results, and write an audit record. |
| `06_balanced_evaluation.R` | Summarize the sensitivity/specificity balanced scheme using the validation-set Youden threshold; output the model comparison, subgroup performance, threshold strategies, and corresponding figures. This step does not retrain the model. |

## Application and Function Checks

| Script | Function |
|---|---|
| `app_predict.R` | Load and verify the saved model and thresholds; process MCV, PDW, and monocyte-percentage inputs; and perform the proportion conversion, model prediction, and classification. |
| `app.R` | Provide the local Shiny interface: accept the three inputs and display the model score, classification tendency, and supplementary zone results. |
| `05_smoke_check.R` | Verify that the app predictions match the saved results, and check the input conversion, threshold classification, and Shiny server-side functionality. |

## Analysis Order

`00_prepare.R` → `01_model_pipeline.R` → `03_supplement.R` (which includes `02_evaluate.R`) → `04_figures.R` → `06_balanced_evaluation.R`.

`05_integrity_check.R` verifies the analysis outputs after `04_figures.R`; `05_smoke_check.R` verifies the application after `06_balanced_evaluation.R`.

## Outputs

Main outputs are written to the working directory:

- `results/` — statistical results and predictions
- `models/` — saved models
- `figures/` — individual figure panels (PDF and PNG)
- `balanced_scheme/` — balanced-scheme results and figures

## Requirements

- R >= 4.0 (analysis performed in R 4.5.1)
- Key packages: `rsample`, `caret`, `recipes`, `themis`, `Boruta`, `glmnet`, `pROC`, `xgboost`, `ranger`, `gbm`, `nnet`, `rpart`, `kernlab`, `ggplot2`, `shiny`
- The source dataset `TB_Pneumonia.csv` placed in the working directory

## Interactive Web App

Research decision-support prototype (not for clinical use):

[https://nana2379723224.shinyapps.io/tb-pneumonia-research/](https://nana2379723224.shinyapps.io/tb-pneumonia-research/)

## Citation

If you use this code in your research, please cite:

> [Author names]. Explainable machine learning triage of tuberculosis versus pneumonia using routine blood parameters. [Journal], [Year]. DOI: [doi]

The code is archived on Zenodo:

> **Concept DOI (all versions)**: [10.5281/zenodo.21151460](https://doi.org/10.5281/zenodo.21151460)

## License

Research and educational purposes only. Not for clinical use.
