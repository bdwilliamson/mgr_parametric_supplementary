# Replicate numerical experiments

This file describes how to reproduce numerical experiments. All analyses were implemented in the freely available R programming language, specifically, R version 4.0.2 or later. 

The files are:
* `00_sim_utils.R`: a function to obtain true parameter values
* `01_gen_data.R`: functions for generating data for simulations
* `02_run_ests.R`: functions to generate a dataset and run all requested estimators for a single simulation replication
* `03_main.R`: run the simulation a set number of times with given arguments
* `04_plots.R`: create plots and tables with results of simulations
* `run_estimators_all_n.bat`: run code on `cmd` for all sample sizes for a given set of estimators