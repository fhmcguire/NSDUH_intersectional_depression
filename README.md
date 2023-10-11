# Depression at the intersection of race/ethnicity, sex/gender, and sexual orientation in a nationally representative sample of US adults: A design-weighted intersectional MAIHDA
This repository contains data and code used to support the manuscript "Depression at the intersection of race/ethnicity, sex/gender, and sexual orientation in a nationally representative sample of US adults: A design-weighted MAIHDA."

**This manuscript is currently undergoing peer review.** A preprint version (Date: 04/13/2023) is available at: https://www.medrxiv.org/content/10.1101/2023.04.13.23288529v1

## Manuscript abstract

This study examined how race/ethnicity, sex/gender, and sexual orientation intersect to socially pattern depression among US adults. We used cross-sectional data from the 2015-2020 National Survey on Drug Use and Health (NSDUH; n=234,722) to conduct design-weighted multilevel analysis of individual heterogeneity and discriminatory accuracy (MAIHDA) for past-year and lifetime major depressive episode (MDE). With 42 intersectional groups constructed from seven race/ethnicity, two sex/gender, and three sexual orientation categories, we estimated age-standardized prevalence and excess/reduced prevalence attributable to intersectional effects (i.e., two-way or higher identity variable interactions). Models revealed heterogeneity between intersectional groups, with prevalence estimates ranging from 1.9–19.7% (past-year) and 4.5–36.5% (lifetime). The intersectional group structure (i.e., general contextual effect) explained 12.8% (past-year) and 12.6% (lifetime) of model variance, indicating key relevance of intersectional groups in determining the population distribution of depression. Main effects indicated, on average, people who were White, women, gay/lesbian, or bisexual had greater odds of MDE. While main effects explained most between-group variance, intersectional effects (past-year: 10.1%; lifetime: 16.5%) indicated heterogeneity around averages, such that groups experienced excess/reduced prevalence compared to main effects expectations. Notably, we extend the MAIHDA framework to calculate nationally representative estimates from complex sample survey data using design-weighted, Bayesian methods.

## Overview of repository structure/contents

1. Original data downloaded from NSDUH are available in Stata format in the "/01_original_data" folder. Files were too large to be uploaded directly, so they are stored in zip folders. This data is processed/cleaned using the "01_nsduh_dataManagement.R" file to create the final data file ("nsduh.RDS") used for analysis.
2. Ideally, MAIHDA model fits would be stored in a "/02_fits" folder. However, these files are too large to be stored on GitHub, so please inquire directly to request access to the model fits.
3. All R code for this project is stored in the main folder and is separated into six files:
    * 00_helper_functions.R (store R functions that automate specific tasks)
    * 01_nsduh_dataManagement.R (data management/cleaning tasks)
    * 02_nsduh_analysis.R (conduct design-weighted MAIHDA [main analysis] and unweighted MAIHDA [senstivity analysis])
    * 03_nsduh_estimates.Rmd (generate and summarize model estimates)
    * 04_nsduh_auc.R (create area under the receiver operating characteristic curve [AUC] plots)
    * 05_tables_figures.Rmd (generate tables/figures for the manuscript)
