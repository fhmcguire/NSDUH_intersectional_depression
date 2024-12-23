# ----------------------------------------------------------------------------
#
# Script name: 04_nsduh_auc.R
#
# Purpose of script: Create area under the receiver operator characteristic
#                    curve (AUC) plots. 
#
#                    Conduct analysis separately for design-weighted and 
#                    unweighted baseline models.
#
# Author: F. Hunter McGuire, MPH
#
# Date created: 2023-09-01
# Last update:  2024-12-23
#
# -----------------------------------------------------------------------------

#############
### SETUP ###
#############

# Load packages
library(pROC) # produce ROC-AUC
library(ggplot2) # Load `tidyverse` suite of packages
library(plyr) # Tools for Splitting, Applying and Combining Data
library(csSampling) # Complex survey sampling with Bayesian methods
library(ggthemes) # Extra themes and functions for `ggplot2`
library(extrafont) # Extra fonts for creating charts/figures
library(rstudioapi) # Retrieve information about an RStudio Edit

## set working directory to source file location
### get full file path
Rpath <- rstudioapi::getSourceEditorContext()$path
### remove file name from file path
Rpath <- gsub('/04_nsduh_auc.R','', Rpath)
### set the working directory
setwd(Rpath)


###############################################################################
########################### Weighted analysis #################################
###############################################################################

# load the helper functions used to extract/manipulate model estimates.
source("00_helper_functions.R")

# load the full dataset
nsduh <- readRDS(file = "nsduh.rds")
# subset the dataset to obs with complete data for analysis 1 (lifetime MDE)
nsduh_ltMDE <- subset(nsduh, sp_1 == 1)
# subset the dataset to obs with complete data for analysis 2 (past-year MDE)
nsduh_pyMDE <- subset(nsduh, sp_2 == 1)

# load the Model 1 fit objects
ltMDE <- readRDS(file = "02_fits/weighted_csSampling/ltMDE_m1.rds")
ltMDE <- ltMDE$stan_fit
pyMDE <- readRDS(file = "02_fits/weighted_csSampling/pyMDE_m1.rds")
pyMDE <- pyMDE$stan_fit

# -----------------------------------------------------------------------------

#########################
### EXTRACT ESTIMATES ###
#########################

# get table of predicted probability (prevalence) estimates for each group

# Define population age distribution weights
age_weights <- c(0.138, 0.159, 0.246, 0.456)
# Note: these percentages come from Table 1
# 18-25 = 75,591 (13.8%)
# 26-34 = 48,242 (15.9%)
# 35-49 = 62,110 (24.6%)
# 50+   = 48,779 (45.6%)

# Generate age-standardized prevalence
ltMDE_pred <-
  age_adj_prevalence(model = ltMDE, 
                     parameters =  c("b_Intercept", "b[1]", "b[2]", "b[3]"),
                     age_weights = age_weights) %>% 
  mutate(pred = as.numeric(Mean) / 100) %>% 
  select(pred, strata = GroupID)

pyMDE_pred <-
  age_adj_prevalence(model = pyMDE, 
                     parameters =  c("b_Intercept", "b[1]", "b[2]", "b[3]"),
                     age_weights = age_weights) %>% 
  mutate(pred = as.numeric(Mean) / 100) %>% 
  select(pred, strata = GroupID)

# merge the predicted probabilities with the main data frames
nsduh_ltMDE_pred <- merge(x = nsduh_ltMDE, 
                          y = ltMDE_pred, 
                          by = "strata", 
                          all.x = TRUE)
nsduh_pyMDE_pred <- merge(x = nsduh_pyMDE, 
                          y = pyMDE_pred, 
                          by = "strata", 
                          all.x = TRUE)

# -----------------------------------------------------------------------------

#####################
### PREPARE PLOTS ###
#####################

# build the ROC curves
pyMDE_roc <- roc(nsduh_pyMDE_pred$pyMDE, 
                 nsduh_pyMDE_pred$pred)
ltMDE_roc <- roc(nsduh_ltMDE_pred$ltMDE,
                 nsduh_ltMDE_pred$pred)

# calculate AUC and 95% confidence interval for each plot
# bootstrap sampling with n = 2000 samples
# round to 3 digits
# past-year MDE
set.seed(828)
pyMDE_auc <- format(round(auc(pyMDE_roc),3), nsmall = 3)
pyMDE_ci <- ci.auc(pyMDE_roc, conf.level = 0.95,
                   method = "bootstrap", boot.n = 2000)
pyMDE_lower <- format(round(pyMDE_ci[1], 3), nsmall = 3)
pyMDE_upper <- format(round(pyMDE_ci[3], 3), nsmall = 3)
# lifetime MDE
set.seed(828)
ltMDE_auc <- format(round(auc(ltMDE_roc),3), nsmall = 3)
ltMDE_ci <- ci.auc(ltMDE_roc, conf.level = 0.95,
                   method = "bootstrap", boot.n = 2000)
ltMDE_lower <- format(round(ltMDE_ci[1], 3), nsmall = 3)
ltMDE_upper <- format(round(ltMDE_ci[3], 3), nsmall = 3)
# store the AUC values in text vector (to be used in the figure)
pyMDE_legend_text <- paste0(
  "AUC = ", pyMDE_auc, " (95% CI: ", pyMDE_lower, ", ", pyMDE_upper, ")")
ltMDE_legend_text <- paste0(
  "AUC = ", ltMDE_auc, " (95% CI: ", ltMDE_lower, ", ", ltMDE_upper, ")")

# -----------------------------------------------------------------------------

############################
### DISPLAY/EXPORT PLOTS ###
############################

# past-year MDE
pROC::ggroc(pyMDE_roc, color = "black", size = 3) +
  ggtitle("Past-year major depressive episode") +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
               color="dimgray", 
               linetype="dashed",
               linewidth = 1.5) +
  labs(x = "Specificity", y = "Sensitivity") +
  ggplot2::annotate("text", x = 0.4, y = 0.1, size = 8,
                    label = pyMDE_legend_text) +
  theme_clean(base_size = 21)
# save the plot
ggsave("03_tables_figures/roc_pyMDE.tiff",
       units=c("in"),
       width = 9,
       height = 7.5,
       dpi = 400)

# lifetime MDE
pROC::ggroc(ltMDE_roc, color = "black", size = 3) +
  ggtitle("Lifetime major depressive episode") +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
               color="dimgray", 
               linetype="dashed",
               linewidth = 1.5) +
  labs(x = "Specificity", y = "Sensitivity") +
  ggplot2::annotate("text", x = 0.4, y = 0.1, size = 8,
                    label = ltMDE_legend_text) +
  theme_clean(base_size = 21)
# save the plot
ggsave("03_tables_figures/roc_ltMDE.tiff",
       units=c("in"),
       width = 9,
       height = 7.5,
       dpi = 400)




###############################################################################
########################### Unwighted analysis ###############################
###############################################################################

# clear the Global Environment
rm(list=ls())

# load the helper functions used to extract/manipulate model estimates.
source("00_helper_functions.R")

# load the full dataset
nsduh <- readRDS(file = "nsduh.rds")
# subset the dataset to obs with complete data for analysis 1 (lifetime MDE)
nsduh_ltMDE <- subset(nsduh, sp_1 == 1)
# subset the dataset to obs with complete data for analysis 2 (past-year MDE)
nsduh_pyMDE <- subset(nsduh, sp_2 == 1)

# load the Model 1 fit objects
ltMDE <- readRDS(file = "02_fits/unweighted_brms/ltMDE_m1_unweighted.rds")
ltMDE <- ltMDE$fit
pyMDE <- readRDS(file = "02_fits/unweighted_brms/pyMDE_m1_unweighted.rds")
pyMDE <- pyMDE$fit

# -----------------------------------------------------------------------------

#########################
### EXTRACT ESTIMATES ###
#########################

# get table of predicted probability (prevalence) estimates for each group
ltMDE_pred = predicted(fit = ltMDE, 
                       n_iter = 8000,
                       int_name = "b_Intercept",
                       group_names = strataName)
pyMDE_pred = predicted(fit = pyMDE, 
                       n_iter = 8000,
                       int_name = "b_Intercept",
                       group_names = strataName)

# merge the predicted probabilities with the main data frames
nsduh_ltMDE_pred <- merge(x = nsduh_ltMDE, 
                          y = ltMDE_pred, 
                          by = "strata", 
                          all.x = TRUE)
nsduh_pyMDE_pred <- merge(x = nsduh_pyMDE, 
                          y = pyMDE_pred, 
                          by = "strata", 
                          all.x = TRUE)

# -----------------------------------------------------------------------------

#####################
### PREPARE PLOTS ###
#####################

# build the ROC curves
pyMDE_roc <- roc(nsduh_pyMDE_pred$pyMDE, 
                 nsduh_pyMDE_pred$predicted)
ltMDE_roc <- roc(nsduh_ltMDE_pred$ltMDE,
                 nsduh_ltMDE_pred$predicted)

# calculate AUC and 95% confidence interval for each plot
# bootstrap sampling with n = 2000 samples
# round to 3 digits
# past-year MDE
set.seed(828)
pyMDE_auc <- format(round(auc(pyMDE_roc),3), nsmall = 3)
pyMDE_ci <- ci.auc(pyMDE_roc, conf.level = 0.95,
                   method = "bootstrap", boot.n = 2000)
pyMDE_lower <- format(round(pyMDE_ci[1], 3), nsmall = 3)
pyMDE_upper <- format(round(pyMDE_ci[3], 3), nsmall = 3)
# lifetime MDE
set.seed(828)
ltMDE_auc <- format(round(auc(ltMDE_roc),3), nsmall = 3)
ltMDE_ci <- ci.auc(ltMDE_roc, conf.level = 0.95,
                   method = "bootstrap", boot.n = 2000)
ltMDE_lower <- format(round(ltMDE_ci[1], 3), nsmall = 3)
ltMDE_upper <- format(round(ltMDE_ci[3], 3), nsmall = 3)
# store the AUC values in text vector (to be used in the figure)
pyMDE_legend_text <- paste0(
  "AUC = ", pyMDE_auc, " (95% CI: ", pyMDE_lower, ", ", pyMDE_upper, ")")
ltMDE_legend_text <- paste0(
  "AUC = ", ltMDE_auc, " (95% CI: ", ltMDE_lower, ", ", ltMDE_upper, ")")

# -----------------------------------------------------------------------------

############################
### DISPLAY/EXPORT PLOTS ###
############################

# past-year MDE
pROC::ggroc(pyMDE_roc, color = "black", size = 3) +
  ggtitle("Past-year major depressive episode") +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
               color="dimgray", 
               linetype="dashed",
               linewidth = 1.5) +
  labs(x = "Specificity", y = "Sensitivity") +
  annotate("text", x = 0.4, y = 0.1, size = 8,
           label = pyMDE_legend_text) +
  theme_clean(base_size = 21)
# save the plot
ggsave("03_tables_figures/unweighted/roc_pyMDE_uw.tiff",
       units=c("in"),
       width = 9,
       height = 7.5,
       dpi = 300)

# lifetime MDE
pROC::ggroc(ltMDE_roc, color = "black", size = 3) +
  ggtitle("Lifetime major depressive episode") +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1),
               color="dimgray", 
               linetype="dashed",
               linewidth = 1.5) +
  labs(x = "Specificity", y = "Sensitivity") +
  annotate("text", x = 0.4, y = 0.1, size = 8,
           label = ltMDE_legend_text) +
  theme_clean(base_size = 21)
# save the plot
ggsave("03_tables_figures/unweighted/roc_ltMDE_uw.tiff",
       units=c("in"),
       width = 9,
       height = 7.5,
       dpi = 300)
