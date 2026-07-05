# TB-Pneumonia-ML-triage

**Machine learning-based triage of tuberculosis versus pneumonia using routine complete blood count parameters.**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21151461.svg)](https://doi.org/10.5281/zenodo.21151461)

## Overview

This repository contains the complete analysis code for an explainable machine learning model that distinguishes pulmonary tuberculosis (TB) from community-acquired pneumonia using five routine complete blood count (CBC) parameters:

- **RBC** — Red Blood Cell Count
- **PDW** — Platelet Distribution Width
- **MONO%** — Monocyte Percentage
- **MPV** — Mean Platelet Volume
- **LMR** — Lymphocyte-to-Monocyte Ratio

## Key Results

| Metric | Value |
|---|---|
| Best AUC (Validation) | 0.866 |
| Sensitivity (Test) | 85.7% |
| Specificity (Test) | 83.3% |
| Model | XGBoost (nrounds = 150, max_depth = 3, eta = 0.05) |
| Triple threshold | Low risk < 14.6%, Gray zone 14.6–26.5%, High risk > 26.5% |

## Repository Structure

```
TB-Pneumonia-ML-triage/
├── R/                          # R analysis scripts
│   ├── 01_data_prep.R          # Data loading, cleaning, Table 1
│   ├── 02_exploratory_figures.R # Figure 2: EDA visualizations
│   ├── 03_feature_selection.R  # Figure 2E, 3A, 3B: Boruta + LASSO + AUC
│   ├── 04_model_benchmarking.R # Figure 3C: 9-model CV comparison
│   ├── 05_final_model.R        # Figure 3D, 3E: XGBoost + ROC + CM
│   ├── 06_calibration_dca.R    # Figure 4A, 4B: Calibration + DCA
│   ├── 07_shap_analysis.R      # Figure 5A, 5B: SHAP beeswarm + force
│   └── 08_shiny_app.R          # Interactive Shiny dashboard
├── data/                       # Raw data (not distributed)
├── docs/
│   └── CODE_AVAILABILITY.md    # Code availability statement
├── .gitignore
└── README.md
```

## Requirements

- R >= 4.0
- Key packages: `caret`, `xgboost`, `pROC`, `glmnet`, `Boruta`, `shapviz`, `dcurves`, `shiny`, `echarts4r`

## Reproduction

1. Clone this repository.
2. Place the raw dataset (contact corresponding author for access) in `data/`.
3. Run scripts sequentially from `01_data_prep.R` to `07_shap_analysis.R`.
4. For the interactive dashboard, place the pre-trained model in the app directory and run `08_shiny_app.R`.

## Interactive Web App

Try the online diagnostic assistant:

[https://0pwfel-0-0.shinyapps.io/TB_AI_App/](https://0pwfel-0-0.shinyapps.io/TB_AI_App/)

## Citation

If you use this code in your research, please cite:

> [Author names], "Explainable machine learning triage of tuberculosis versus pneumonia using routine blood parameters," [Journal], [Year]. DOI: [doi]

## License

Research and educational purposes only. Not for clinical use.
