## ----------------------------------------------------------------------------
##
## Script name:       02_nsduh_analysis.R
##
## Purpose:           Conduct design-weighted intersectional MAIHDA for
##                    past-year and lifetime major depressive episode (MDE)
##                    using data from NSDUH 2015-2020. Intersectional groups
##                    are defined by race/ethnicity, sex/gender, and sexual
##                    orientation. All models are adjusted for age category. 
##                    Models are fit with Bayesian two-level logisitic models
##                    that account for NSDUH complex sampling design weights.
##
##                    Conduct sensitivity analysis with unweighted Bayesian
##                    two-level logistic models fit with brms.
##
## Author:            F. Hunter McGuire, PhD, MPH
##
## Date created:      2023-08-17
## Last update:       2023-10-11
##
## ----------------------------------------------------------------------------
##
## Data source: National Survey on Drug Use and Health (NSDUH) 2015-2020
##
## Input data: nsduh.rds
##
## Output data: Model fits (pyMDE_m1, pyMDE_m2, ltMDE_m1, ltMDE_m2)
##
## Variables:
##
##    'ltMDE'      = lifetime major depressive episode (0 = no, 1 = yes)
##    'pyMDE'      = past-year major depressive episode (0 = no, 1 = yes)
##    'mhservice'  = past-year mental health service use (0 = no, 1 = yes)
##
##    'asian'      = race indicator (1 = Asian)
##    'black'      = race indicator (1 = Black / African-American)
##    'hispanic'   = race indicator (1 = Hispanic)
##    'naan'       = race indicator (1 = Native American / Alaskan Native)
##    'nhpi'       = race indicator (1 = Native Hawaiian / Pacific Islander)
##    'multi'      = race indicator (1 = Multiracial)
##    'white'      = race indicator (1 = White)
##
##    'woman'      = gender indicator (1 = Woman)
##    'man'        = gender indicator (1 = Man)
##
##    'hetero'     = sexual orientation indicator (1 = Heterosexual)
##    'gay'        = sexual orientation indicator (1 = Gay / Lesbian)
##    'bisexual'   = sexual orientation indicator (1 = Bisexual)
##
##    'age1217'    = aged 12-17 years old (0 = no, 1 = yes)
##    'age1825'    = aged 18-25 years old (0 = no, 1 = yes)
##    'age2634'    = aged 26-34 years old (0 = no, 1 = yes)
##    'age3549'    = aged 35-49 years old (0 = no, 1 = yes)
##    'age50plus'  = aged 50+ years old (0 = no, 1 = yes)
##    
##    'sp_1'       = subpop for lifetime MDE analysis (0 = no, 1 = yes)
##    'sp_2'       = subpop for past-year MDE analysis (0 = no, 1 = yes)
##
##    'verep'      = NSDUH cluster weights
##    'vestr'      = NSDUH strata weights
##    'adjwt6'     = ADJUSTED NSDUH person-level weights for number of data 
##                   collection years (6 years: 2015-2020)
##    'wt_norm'    = scaled person-level weights (mean=1); input='adjwt6'
##
## ----------------------------------------------------------------------------


# set working directory to source file location
### get full file path
RMDpath = rstudioapi::getSourceEditorContext()$path
### remove file name from file path
RMDpath = gsub('/02_nsduh_analysis_update.R','', RMDpath)
### set the working directory
setwd(RMDpath)

# Load packages
library(csSampling)
library(dplyr)
library(brms)
library(rstan)
library(survey)
library(cmdstanr)

# Load data (cleaned in the "nsduh_dataManagement.R" file)
nsduh <- readRDS(file="nsduh.RDS")

# Define priors
nsduh_priors <- c(
  # logit scale
  prior(normal(0, 2), class = Intercept),
  prior(normal(0, 2), class = b),
  # variance components
  prior(inv_gamma(2, 2), class = sd))



####################################
#### Setup survey design object ####
####################################

# Scale the sampling weights so the mean=1
nsduh$wt_norm <- nsduh$adjwt6*(length(nsduh$adjwt6)/sum(nsduh$adjwt6)) 

# Create overall survey design object
design <- 
  svydesign(
    id = ~ verep, # cluster weight
    strata = ~ vestr, # strata weight
    weights = ~ wt_norm, # scaled person-level weight
    nest = TRUE,
    data = nsduh
  )

# Create replicate design object
repdesign <- survey::as.svrepdesign(design = design,
                                    type = "mrbbootstrap",
                                    replicates = 100)


# Create subpopulation (analytic sample) survey design object
# sp_1 = lifetime MDE subpopulation
design_sp1 <- subset(repdesign, sp_1==1) 
# Check mean of weights from subpop survey design object
summary(design_sp1$pweights)
# Calibrate weights so the mean=1
design_sp1 <- calibrate(design_sp1,
                        ~ 1,
                        c(`(Intercept)` = length(design_sp1$pweights)),
                        compress = FALSE)
# Check mean after calibration
summary(design_sp1$pweights) # Mean=1


# sp_2 = past-year MDE subpopulation
design_sp2 <- subset(repdesign, sp_2==1) 
# Check mean of weights from subpop survey design object
summary(design_sp2$pweights)
# Calibrate weights so the mean=1
design_sp2 <- calibrate(design_sp2,
                        ~ 1,
                        c(`(Intercept)` = length(design_sp2$pweights)),
                        compress = FALSE)
# Check mean after calibration
summary(design_sp2$pweights) # Mean=1




#################
# PAST-YEAR MDE #
#################

#################
#### Model 1 ####
#################

# Write and save Stan file
stancode_pyMDE_m1 <- make_stancode(
  brmsformula(
    as.factor(pyMDE)|weights(wt_norm) ~ 
      1 + age2634 + age3549 + age50plus + (1 | strata)
    ), 
  data = nsduh,
  prior = nsduh_priors,
  family = bernoulli(link="logit"),
  save_model = "02_fits/weighted_stancode/brms_pyMDE_m1.stan")
# Load the Stan file
modbrms_pyMDE_m1 <- stan_model("02_fits/weighted_stancode/brms_pyMDE_m1.stan")
# Set up data
databrms_pyMDE_m1 <- make_standata(
  brmsformula(
    as.factor(pyMDE)|weights(wt_norm) ~ 
      1 + age2634 + age3549 + age50plus + (1 | strata)
    ),
  data = design_sp2$variables,
  prior = nsduh_priors,
  family = bernoulli(link="logit"))
# Set Stan model weights = survey design weights
databrms_pyMDE_m1$weights <- design_sp2$pweights

# Model estimation
set.seed(828) # for reproducible results
pyMDE_m1 <- cs_sampling(
  svydes = design_sp2, # use the subpop-defined survey object
  mod_stan = modbrms_pyMDE_m1,
  data_stan = databrms_pyMDE_m1,
  ctrl_stan = list(chains = 4,
                   iter = 4000, # 16,000 total iterations (8,000 post-warmup)
                   warmup = 2000, # 8,000 warmups
                   prior = nsduh_priors,
                   backend = "cmdstanr",
                   threads = threading(2),
                   thin = 1),
  rep_design = TRUE,
  sampling_args = list(cores = 4))

# Save the model
saveRDS(pyMDE_m1, file="02_fits/weighted_csSampling/pyMDE_m1.RDS")


#################
#### Model 2 ####
#################

# Write and save Stan file
stancode_pyMDE_m2 <- make_stancode(
  brmsformula(
    as.factor(pyMDE)|weights(wt_norm) ~ 
      1 + age2634 + age3549 + age50plus +
      asian + black + hispanic + naan + nhpi + multi +
      woman + gay + bisexual + (1 | strata)
    ),
  data = nsduh,
  prior = nsduh_priors,
  family = bernoulli(link="logit"),
  save_model = "02_fits/weighted_stancode/brms_pyMDE_m2.stan")
# Load the Stan file
modbrms_pyMDE_m2 <- stan_model("02_fits/weighted_stancode/brms_pyMDE_m2.stan")
# Set up data
databrms_pyMDE_m2 <- make_standata(
  brmsformula(
    as.factor(pyMDE)|weights(wt_norm) ~ 
      1 + age2634 + age3549 + age50plus +
      asian + black + hispanic + naan + nhpi + multi +
      woman + gay + bisexual + (1 | strata)
    ),
  data = design_sp2$variables,
  prior = nsduh_priors,
  family = bernoulli(link="logit"))
# Set Stan model weights = survey design weights
databrms_pyMDE_m2$weights <- design_sp2$pweights


# Model estimation
set.seed(828) # for reproducible results
pyMDE_m2 <- cs_sampling(
  svydes = design_sp2, # use the subpop-defined survey object
  mod_stan = modbrms_pyMDE_m2,
  data_stan = databrms_pyMDE_m2,
  ctrl_stan = list(chains = 4,
                   iter = 4000, # 16,000 total iterations (8,000 post-warmup)
                   warmup = 2000, # 8,000 warmups
                   prior = nsduh_priors,
                   backend = "cmdstanr",
                   threads = threading(2),
                   thin = 1),
  rep_design = TRUE,
  sampling_args = list(cores = 8))

# Save the model
saveRDS(pyMDE_m2, file="02_fits/weighted_csSampling/pyMDE_m2.RDS")




#################
# LIFETIME MDE ##
#################

#################
#### Model 1 ####
#################

# Write and save Stan file
stancode_ltMDE_m1 <- make_stancode(
  brmsformula(
    as.factor(ltMDE)|weights(wt_norm) ~ 
      1 + age2634 + age3549 + age50plus + (1 | strata)
  ), 
  data = nsduh,
  prior = nsduh_priors,
  family = bernoulli(link="logit"),
  save_model = "02_fits/weighted_stancode/brms_ltMDE_m1.stan")
# Load the Stan file
modbrms_ltMDE_m1 <- stan_model("02_fits/weighted_stancode/brms_ltMDE_m1.stan")
# Set up data
databrms_ltMDE_m1 <- make_standata(
  brmsformula(
    as.factor(ltMDE)|weights(wt_norm) ~ 
      1 + age2634 + age3549 + age50plus + (1 | strata)
  ),
  data = design_sp1$variables,
  prior = nsduh_priors,
  family = bernoulli(link="logit"))
# Set Stan model weights = survey design weights
databrms_ltMDE_m1$weights <- design_sp1$pweights

# Model estimation
set.seed(828) # for reproducible results
ltMDE_m1 <- cs_sampling(
  svydes = design_sp1, # use the subpop-defined survey object
  mod_stan = modbrms_ltMDE_m1,
  data_stan = databrms_ltMDE_m1,
  ctrl_stan = list(chains = 4,
                   iter = 4000, # 16,000 total iterations (8,000 post-warmup)
                   warmup = 2000, # 8,000 warmups
                   prior = nsduh_priors,
                   backend = "cmdstanr",
                   threads = threading(2),
                   thin = 1),
  rep_design = TRUE,
  sampling_args = list(cores = 8))

# Save the model
saveRDS(ltMDE_m1, file="02_fits/weighted_csSampling/ltMDE_m1.RDS")


#################
#### Model 2 ####
#################

# Write and save Stan file
stancode_ltMDE_m2 <- make_stancode(
  brmsformula(
    as.factor(ltMDE)|weights(wt_norm) ~ 
      1 + age2634 + age3549 + age50plus +
      asian + black + hispanic + naan + nhpi + multi +
      woman + gay + bisexual + (1 | strata)
  ),
  data = nsduh,
  prior = nsduh_priors,
  family = bernoulli(link="logit"),
  save_model = "02_fits/weighted_stancode/brms_ltMDE_m2.stan")
# Load the Stan file
modbrms_ltMDE_m2 <- stan_model("02_fits/weighted_stancode/brms_ltMDE_m2.stan")
# Set up data
databrms_ltMDE_m2 <- make_standata(
  brmsformula(
    as.factor(ltMDE)|weights(wt_norm) ~ 
      1 + age2634 + age3549 + age50plus +
      asian + black + hispanic + naan + nhpi + multi +
      woman + gay + bisexual + (1 | strata)
  ),
  data = design_sp1$variables,
  prior = nsduh_priors,
  family = bernoulli(link="logit"))
# Set Stan model weights = survey design weights
databrms_ltMDE_m2$weights <- design_sp1$pweights


# Model estimation
set.seed(828) # for reproducible results
ltMDE_m2 <- cs_sampling(
  svydes = design_sp1, # use the subpop-defined survey object
  mod_stan = modbrms_ltMDE_m2,
  data_stan = databrms_ltMDE_m2,
  ctrl_stan = list(chains = 4,
                   iter = 4000, # 16,000 total iterations (8,000 post-warmup)
                   warmup = 2000, # 8,000 warmups
                   prior = nsduh_priors,
                   backend = "cmdstanr",
                   threads = threading(2),
                   thin = 1),
  rep_design = TRUE,
  sampling_args = list(cores = 8))

# Save the model
saveRDS(ltMDE_m2, file="02_fits/weighted_csSampling/ltMDE_m2.RDS")




###############################################################################
##################### Sensitivity: Unweighted analysis ########################
###############################################################################

# remove all objects from global environment
rm(list=ls())

# Load packages
library(brms)
library(cmdstanr)
library(dplyr)
library(marginaleffects)

# Load data (cleaned in the "nsduh_dataManagement.R" file)
nsduh <- readRDS(file="nsduh.RDS")

# subset data to analytic sample
nsduh_pyMDE <- subset(nsduh, sp_2==1)
nsduh_ltMDE <- subset(nsduh, sp_1==1)

# Define priors
nsduh_priors <- c(
  # logit scale
  prior(normal(0, 2), class = Intercept),
  prior(normal(0, 2), class = b),
  # variance components
  prior(inv_gamma(2, 2), class = sd))


#####################
### PAST-YEAR MDE ###
#####################

# age-adjusted ("null") model
pyMDE_m1_unweighted <- 
  brm(pyMDE ~ 1 + age2634 + age3549 + age50plus + (1 | strata),
      data = nsduh_pyMDE,
      family = bernoulli(link="logit"),
      prior = nsduh_priors,
      chains = 4,
      iter = 4000,
      warmup = 2000,
      thin = 1,
      threads = threading(2),
      cores = 4,
      backend = "cmdstanr")
# save the model
saveRDS(pyMDE_m1_unweighted, 
        "02_fits/unweighted_brms/pyMDE_m1_unweighted.rds")


# fully-adjusted ("intersectional interaction") model
pyMDE_m2_unweighted <- 
  brm(pyMDE ~ 1 + age2634 + age3549 + age50plus + 
        asian + black + hispanic + naan + nhpi + multi +
        woman + gay + bisexual + (1 | strata),
      data = nsduh_pyMDE,
      family = bernoulli(link="logit"),
      prior = nsduh_priors,
      chains = 4,
      iter = 4000,
      warmup = 2000,
      thin = 1,
      threads = threading(2),
      cores = 4,
      backend = "cmdstanr")
# save the model
saveRDS(pyMDE_m2_unweighted, 
        "02_fits/unweighted_brms/pyMDE_m2_unweighted.rds")


####################
### LIFETIME MDE ###
####################

# age-adjusted ("null") model
ltMDE_m1_unweighted <- 
  brm(ltMDE ~ 1 + age2634 + age3549 + age50plus + (1 | strata),
      data = nsduh_ltMDE,
      family = bernoulli(link="logit"),
      prior = nsduh_priors,
      chains = 4,
      iter = 4000,
      warmup = 2000,
      thin = 1,
      threads = threading(2),
      cores = 4,
      backend = "cmdstanr")
# save the model
saveRDS(ltMDE_m1_unweighted, 
        "02_fits/unweighted_brms/ltMDE_m1_unweighted.rds")


# fully-adjusted ("intersectional interaction") model
ltMDE_m2_unweighted <- 
  brm(ltMDE ~ 1 + age2634 + age3549 + age50plus + 
        asian + black + hispanic + naan + nhpi + multi +
        woman + gay + bisexual + (1 | strata),
      data = nsduh_ltMDE,
      family = bernoulli(link="logit"),
      prior = nsduh_priors,
      chains = 4,
      iter = 4000,
      warmup = 2000,
      thin = 1,
      threads = threading(2),
      cores = 4,
      backend = "cmdstanr")
# save the model
saveRDS(ltMDE_m2_unweighted, 
        "02_fits/unweighted_brms/ltMDE_m2_unweighted.rds")

