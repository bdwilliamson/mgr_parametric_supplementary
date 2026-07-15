# Replicate data analyses

This file describes how to reproduce the data analysis. The data are not provided to protect confidentiality. All analyses were implemented in the freely available R programming language, specifically, R version 4.0.2 or later. 

The files are:
* `00_vccc_cleaning_utils.R`: utility functions for cleaning the VCCC data
* `01_clean_vccc_data.R`: clean the VCCC data by combining correct datasets together, creating derived variables (including binary outcome and exposure)
* `02_sample_phase_2.R`: sample 1000 two-phase datasets using case-control sampling
* `03_main.R`: run all estimators, including ideal-data estimators
* `04_plots_and_tables.R`: create plots and tables with results