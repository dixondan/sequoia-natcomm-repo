# Sequoia Mortality Analysis

This repository contains the workflow used to estimate giant sequoia mortality and evaluate prescribed fire scenarios.

## Run order

Rscript r/models/0_calibrate-gamma.R  
Rscript r/models/2_simulate-counts.R  
Rscript r/models/3_model-main.R 500 6 FALSE  
Rscript r/models/5_scenarios.R  

Then in R:
source("r/models/4_plot-effects.R")

## Inputs

data/intermediate/model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv  
data/intermediate/stage2-bern_matrix_gamma_6000.rds  

## Outputs

results/tables/  
results/figures/  

## Notes

Python scripts were used for preprocessing but are not required to run the analysis.