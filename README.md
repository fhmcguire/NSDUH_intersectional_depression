# Depression at the intersection of race/ethnicity, sex/gender, and sexual orientation in a nationally representative sample of US adults: A design-weighted intersectional MAIHDA
This repository contains data and code used to support the manuscript "Depression at the intersection of race/ethnicity, sex/gender, and sexual orientation in a nationally representative sample of US adults: A design-weighted MAIHDA."

As of June 14, 2024, this manuscript has been published in the **[American Journal of Epidemiology](https://academic.oup.com/aje/advance-article/doi/10.1093/aje/kwae121/7693604?utm_source=authortollfreelink&utm_campaign=aje&utm_medium=email&guestAccessKey=6aeb2470-7617-41eb-9529-5ab0bbfb3804)**. 

## Manuscript abstract

This study examined how race/ethnicity, sex/gender, and sexual orientation intersect under interlocking systems of oppression to socially pattern depression among US adults. With cross-sectional data from the 2015-2020 National Survey on Drug Use and Health (NSDUH; n=234,722), we conducted design-weighted multilevel analysis of individual heterogeneity and discriminatory accuracy (MAIHDA) under an intersectional framework to predict past-year and lifetime major depressive episode (MDE). With 42 intersectional groups constructed from seven race/ethnicity, two sex/gender, and three sexual orientation categories, we estimated age-standardized prevalence and excess/reduced prevalence attributable to two-way or higher interaction effects. Models revealed heterogeneity across groups, with prevalence ranging from 1.9–19.7% (past-year) and 4.5–36.5% (lifetime). Approximately 12.7% (past-year) and 12.5% (lifetime) of total individual variance were attributable to between-group differences, indicating key relevance of intersectional groups in describing the population distribution of depression. Main effects indicated, on average, people who were White, women, gay/lesbian, or bisexual had greater odds of MDE. Main effects explained most between-group variance. Interaction effects (past-year: 10.1%; lifetime: 16.5%) indicated a further source of heterogeneity around averages with groups experiencing excess/reduced prevalence compared to main effects expectations. We extend the MAIHDA framework to calculate nationally representative estimates from complex sample survey data using design-weighted, Bayesian methods.

## Overview of repository structure/contents

1. Original data downloaded from NSDUH are available in Stata format in the "/01_original_data" folder. Files were too large to be uploaded directly, so they are stored in zip folders. This data is processed/cleaned using the "01_nsduh_dataManagement.R" file to create the final data file ("nsduh.RDS") used for analysis.
2. Ideally, MAIHDA model fits would be stored in a "/02_fits" folder. However, these files are too large to be stored on GitHub, so please inquire directly to request access to the model fits.
3. All R code for this project is stored in the main folder and is separated into six files:
    * 00_helper_functions.R (store R functions that automate specific tasks)
    * 01_nsduh_dataManagement.R (data management/cleaning tasks)
    * 02_nsduh_analysis.R (conduct design-weighted MAIHDA [main analysis] and unweighted MAIHDA [sensitivity analysis])
    * 03_nsduh_estimates.Rmd (generate and summarize model estimates)
    * 04_nsduh_auc.R (create area under the receiver operating characteristic curve [AUC] plots)
    * 05_tables_figures.Rmd (generate tables/figures for the manuscript)
