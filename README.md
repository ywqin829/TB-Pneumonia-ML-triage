# TB-Pneumonia-ML-triage

**Machine learning-based triage of tuberculosis versus pneumonia using routine complete blood count parameters.**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21151461.svg)](https://doi.org/10.5281/zenodo.21151461)

## Overview

This repository contains the analysis code for an explainable machine learning model that distinguishes pulmonary tuberculosis (TB) from pneumonia using three routine complete blood count (CBC) parameters:

- **MCV** — Mean Corpuscular Volume
- **PDW** — Platelet Distribution Width
- **MONO%** — Monocyte Fraction

Feature engineering used Boruta and LASSO regression followed by direction-aware univariate AUC screening and exhaustive feature-combination search; the resulting three-feature panel was used to benchmark nine ML algorithms. The selected eXtreme Gradient Boosting (XGBoost) model was interpreted with SHAP and Decision Curve Analysis and implemented as a local research dashboard.

## Cohort

Retrospective cohort of **2,382 records** (711 confirmed active pulmonary TB; 1,671 pneumonia). Records were partitioned by diagnosis with `rsample` (`seed = 123`) into training (n = 1,428), validation (n = 476), and internal hold-out test (n = 478) sets at an approximate 6:2:2 ratio.

## Key Results

| Metric | Test set |
|---|---|
| AUC | 0.833 (95% CI, 0.793–0.873) |
| Model | XGBoost |
| Panel | MCV, PDW, MONO% |
| Primary cutoff (validation Youden) | 0.422576 |
| Sensitivity | 84.6% (95% CI, 77.8%–89.6%) |
| Specificity | 71.0% (95% CI, 66.0%–75.6%) |
| Balanced accuracy | 77.8% (95% CI, 73.9%–81.5%) |
| Brier score | 0.166 |
| Calibration intercept | −0.825 |

**Operating strategies.** The primary balanced operating point is the validation-derived Youden cutoff (0.422576), adopted after test performance had been reviewed (a disclosed post-result strategy adjustment). A supplementary dual-threshold triage strategy is retained: sensitivity ≥ 90% cutoff at 0.290 (test sensitivity 91.6%, specificity 50.7%) and rule-in specificity ≥ 90% cutoff at 0.634, defining lower (< 0.290), intermediate (0.290–0.634), and higher (≥ 0.634) score zones.

## Repository Structure

```
TB-Pneumonia-ML-triage/
├── R/
│   ├── 00_prepare.R                 # Source audit and frozen 60/20/20 stratified split (seed 123)
│   ├── 01_model_pipeline.R          # Boruta -> LASSO -> AUC -> combinations -> 9-algorithm benchmarking
│   ├── 02_evaluate.R                # Locked validation threshold and one-pass test evaluation
│   ├── 03_supplement.R              # Supplemental descriptive, stability, demographic outputs
│   ├── 04_figures.R                 # Scientific figures rebuilt from locked outputs (per-panel export)
│   ├── 05_integrity_check.R         # Independent, read-only audit of saved outputs
│   ├── 05_smoke_check.R             # Delivery smoke check (predictions, unit conversion, Shiny server)
│   ├── 06_balanced_evaluation.R     # Author-selected balanced operating strategy
│   ├── app.R                        # Shiny research dashboard
│   └── app_predict.R                # Prediction utilities for the dashboard
├── data/                            # Raw data (not distributed)
├── docs/
│   └── CODE_AVAILABILITY.md         # Code availability statement
├── .gitignore
└── README.md
```

## Requirements

- R >= 4.0 (analysis performed in R 4.5.1)
- Key packages: `rsample`, `caret`, `recipes`, `themis`, `Boruta`, `glmnet`, `pROC`, `xgboost`, `ggplot2`, `shiny`

## Reproduction

1. Clone this repository.
2. Place the source dataset as `TB_Pneumonia.csv` in the working directory (contact the corresponding author for access; raw data are not distributed).
3. Run the analysis pipeline in order from the working directory:

   ```
   Rscript R/00_prepare.R
   Rscript R/01_model_pipeline.R
   Rscript R/02_evaluate.R
   Rscript R/03_supplement.R
   Rscript R/04_figures.R
   Rscript R/05_integrity_check.R
   Rscript R/06_balanced_evaluation.R
   ```

   Intermediate outputs are written to `results/`, `figures/`, and `models/`.
4. For the interactive dashboard, run the Shiny app from a directory containing the saved `results/` and `models/` outputs (see `R/app.R`).

## Interactive Web App

Research decision-support prototype (not for clinical use):

[https://nana2379723224.shinyapps.io/tb-pneumonia-research/](https://nana2379723224.shinyapps.io/tb-pneumonia-research/)

## Citation

If you use this code in your research, please cite:

> [Author names], "Explainable machine learning triage of tuberculosis versus pneumonia using routine blood parameters," [Journal], [Year]. DOI: [doi]

## License

Research and educational purposes only. Not for clinical use.
