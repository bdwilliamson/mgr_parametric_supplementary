# Supplementary materials for the parametric MGR paper

This repository contains code to reproduce the analyses in ["Marginal generalized raking with parametric working models"](https://arxiv.org/abs/) by Williamson, Zou, Lumley, Shepherd, and Shaw [(_arXiv_, 2026+)](https://doi.org/10.48550/arXiv.2605.31567). All analyses were implemented in the freely available R programming language, specifically, R version 4.0.2 or later. 

Further, all analyses were performed on multi-core virtual machines running the Windows operating system. There is code provided in each of the directories that will run the simulations from the Windows command prompt (`cmd`). If you use another system or batch job scheduler (e.g., a Linux cluster), these can serve as templates. You may also run the analyses locally, but they will take more time.

This README file provides an overview of the code available in the repository. We have separated our code into two directories based on the two main objectives of the manuscript. Each directory contains a separate README file that describes the code in that directory.

An important file is `00_utils.R`. This file implements all methods.

## The `sims` directory

This directory contains code to reproduce numerical experiments that describe the operating characteristics of our proposed method under varying data-generating mechanisms, in support of our main motivating data analysis. 

## The `data_analysis` directory

This directory contains code to reproduce the main data analysis, estimating the average treatment effect and relative risk of a person living with HIV (PWH) experiencing an AIDS-defining event or death within 5 years of treatment initiation, comparing PWH who initiated a combination antiretroviral therapy (ART) regimen containing a protease inhibitor (PI) against those who initiated a combination ART regiment that did not contain a PI. The data are not provided.

## Issues

If you encounter any bugs or have any specific questions about the analysis, please [file an issue](https://github.com/bdwilliamson/mgr_parametric_supplementary/issues).