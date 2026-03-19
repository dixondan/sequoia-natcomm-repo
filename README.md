# Sequoia Mortality Analysis

This repository contains the workflow used to estimate giant sequoia mortality and evaluate prescribed fire scenarios using a Bayesian tree-level model.

---

## Overview

The analysis combines calibrated mortality probabilities with a hierarchical logistic regression (GLMM) to estimate:

- Observed mortality  
- Mortality without prescribed burns  
- Mortality under universal prescribed burning  

The model operates at the individual tree level with grove-level random effects.

---

## Pipeline Overview

Run scripts in this order from the repository root:

1. Calibration  
2. Mortality simulation  
3. Hierarchical model (ensemble GLMM)  
4. Plotting
5. Scenario analysis  

---

## Run Instructions

### 1. Calibration
Rscript r/models/0_calibrate-gamma.R

### 2. Simulate mortality counts
Rscript r/models/2_simulate-counts.R

### 3. Fit ensemble model
Rscript r/models/3_model-main.R 500 6 FALSE

Arguments:
- 500 = number of ensemble runs  
- 6 = number of CPU cores  
- FALSE = full run (not test mode)  

### 4. Scenario analysis
Rscript r/models/5_scenarios.R

### 5. Plot results (run in R)
source("r/models/4_plot-effects.R")

---

## Required Inputs

The pipeline starts from processed data:

data/intermediate/model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv  
data/intermediate/stage2-bern_matrix_gamma_6000.rds  

---

## Outputs

Tables are written to:
results/tables/

Figures are written to:
results/figures/

---

## Key Results (full run)

Typical results from the full ensemble:

- Observed mortality: ~8,500 trees  
- No prescribed burns: ~10,400 trees  
- Universal prescribed burns: ~4,500 trees  

Interpretation:

- Existing prescribed burns prevented ~1,900 deaths  
- Universal treatment could have prevented ~4,000 deaths  

---

## Model Summary

A Bayesian hierarchical logistic regression is used:

- Response: tree mortality (Bernoulli)  
- Predictors: structure, fuels, topography, prescribed fire  
- Random effects: grove-level variation  
- Ensemble approach: repeated subsampling + posterior aggregation  

---

## Data Availability

Data Availability
Sequoia Tree Inventory (STI)
The tree inventory data (data/raw/AllTreeRecentGrowthAssessments_20230802-KNP_CASTLE.shp) is provided in this repository with permission from Sequoia and Kings Canyon National Parks. This dataset remains the property of the National Park Service. Researchers wishing to use this data independently should contact:
Sequoia and Kings Canyon National Parks
47050 Generals Highway, Three Rivers, CA 93271
https://www.nps.gov/seki/
All other data in this repository is released under the MIT License (see LICENSE).


Data are not included in this repository due to size constraints.

The analysis can be reproduced starting from:

data/intermediate/model-trees-*.csv  
data/intermediate/stage2-bern_matrix_gamma_6000.rds  

These can be archived separately or provided upon request.

--- 

## Runtime Notes

- Full model run (n_iter = 500) may take several hours  
- Small runs (n_iter < 10) are for testing only and will produce unstable results  

---

## Python Preprocessing (coming soon)

Python scripts used for remote sensing preprocessing and feature extraction will be added in a future update.

The current repository reproduces all modeling and scenario results starting from the processed dataset.

---

## Notes

- Run all scripts from the repository root (or use an .Rproj)  
- Do not modify files in data/raw/  
- Intermediate files can be regenerated  
- Final results are written to results/