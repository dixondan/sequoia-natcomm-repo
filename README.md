# Sequoia Mortality Analysis

This repository contains the workflow used to estimate giant sequoia mortality and evaluate survival following prescribed fire scenarios using a Bayesian tree-level model.

---

## Overview

The analysis combines calibrated mortality probabilities with a hierarchical logistic regression (GLMM) framework to estimate:

- Observed mortality
- Mortality without prescribed burns
- Mortality under universal prescribed burning

The model operates at the individual tree level with grove-level random effects.

---

## Repository Structure

```text
python/     Remote sensing preprocessing and feature compilation
r/models/   Calibration, simulation, Bayesian models, and scenarios
results/    Output figures and tables
data/       Input and intermediate data (archived separately)
```

---

## Pipeline Overview

Run scripts in this order from the repository root.

### Python preprocessing

```bash
python python/01_prepare_crowns.py
python python/02_prepare_labels.py
python python/03_pb_upslope.py
python python/03a_lidar.py
python python/03b_topo.py
python python/04_compile_model_data.py
```

### R modeling workflow

#### 0. Calibration

```bash
Rscript r/models/0_calibrate-gamma.R
```

#### 1. Validation

```bash
Rscript r/models/1_cv-gamma.R
```

#### 2. Mortality simulation

```bash
Rscript r/models/2_simulate-counts.R
```

#### 3. Ensemble GLMMs

```bash
Rscript r/models/3_model-main.R 500 6 FALSE
```

Arguments:

- `500` = number of ensemble iterations
- `6` = number of CPU cores
- `FALSE` = full run (not test mode)

#### 4. Scenario analysis

```bash
Rscript r/models/5_scenarios.R
```

#### 5. Plotting

Run in an R session:

```r
source("r/models/4_plot-effects.R")
```

---

## Required Inputs

The modeling workflow starts from:

```text
data/intermediate/model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv
data/intermediate/stage2-bern_matrix_gamma_6000.rds
```

---

## Outputs

Tables are written to:

```text
results/tables/
```

Figures are written to:

```text
results/figures/
```

---

## Model Summary

A Bayesian hierarchical logistic regression is used:

- Response: tree mortality (Bernoulli)
- Predictors: structure, fuels, topography, and prescribed fire history
- Random effects: grove-level variation
- Ensemble framework: repeated spatial subsampling and posterior aggregation

---

## Data Availability

Data are archived separately from this repository.

The modeling workflow can be reproduced starting from:

```text
data/intermediate/model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv
data/intermediate/stage2-bern_matrix_gamma_6000.rds
```

Some preprocessing scripts require lidar, topographic, and fire-history datasets that are not redistributed within this repository.

---

## Runtime Notes

- Full ensemble runs (`n_iter = 500`) may require several hours
- Small runs (`n_iter < 10`) are intended for testing only

---

## Notes

- Run all scripts from the repository root
- Do not modify files in `data/raw/`
- Intermediate files can be regenerated
- Final outputs are written to `results/`