# -----------------------------------------------------------------------------
#
# Script name: 00_helper_functions.R
#
# Author: F. Hunter McGuire, MPH
#
# Date created: 2022-08-27
# Last update:  2023-10-11
#
# -----------------------------------------------------------------------------
#
# This code defines functions and other helper tools used to extract/transform 
# model estimates.
#
# First, we create the `groupName` object to store the intersectional group
# names. We also create the `groupNameShort` object to store a shortened
# version of the names (used in plotting figures).
#
# To examine model convergence diagnostics for each parameter, the 
# `autocor_plots` produces autocorrelation plots and the `trace_plots` 
# function produces trace plots. 
#
# Since logistic regression was used, we used the `logit2prob` function 
# to convert model estimates from the logit to the probability scale. 
#
# The `renameArray` function renames parameter values for the intersectional 
# groups to make the convergence diagnostic plots easier to read/interpret.
#
# For each outcome, the `age_adj_prevalence` calculates the model-predicted,
# age-standardized prevalence of the outcome. It summarizes them as means +
# 95% credible intervals (CI).
#
# For each model, the `variance` function calculates standard deviation (SD), 
# variance, variance partition coefficient (VPC), and proportional change 
# in variance (PCV). As inputs, it takes a list of stanfit (rstan) objects.
#
# For each intersectional group and using the full model (i.e., intersectional
# interaction model), the `calculate_effects` function was used to 
# estimate main effects (i.e., the sum of identity variable coefficients 
# relevant to each group), the intersectional effects (i.e., two-way or higher
# interactions between identity variables), and the full effects (i.e., the 
# main effects and intersectional effects). 
#
# The `add_main_effects` function was used to add the main effects associated
# with each intersectional group (e.g., Black gay woman = Black + Gay + Woman).
# This was specifically used within the `calculate_effects` function.
#
# We used the `summarize_effects` function to calculate means and 
# 95% credible intervals for the total, main, and intersectional effects 
# for each intersectional group. This was specifically used within the 
# `calculate_effects` function.
#
# -----------------------------------------------------------------------------


# Define vector of intersectional group/strata names
groupName <- c("White heterosexual men",
               "White gay men",
               "White bisexual men",
               "White heterosexual women",
               "White gay/lesbian women",
               "White bisexual women",
               "Black heterosexual men",
               "Black gay men",
               "Black bisexual men",
               "Black heterosexual women",
               "Black gay/lesbian women",
               "Black bisexual women",
               "Hispanic/Latine heterosexual men",
               "Hispanic/Latine gay men",
               "Hispanic/Latine bisexual men",
               "Hispanic/Latine heterosexual women",
               "Hispanic/Latine gay/lesbian women",
               "Hispanic/Latine bisexual women",
               "Asian heterosexual men",
               "Asian gay men",
               "Asian bisexual men",
               "Asian heterosexual women",
               "Asian gay/lesbian women",
               "Asian bisexual women",
               "NAAN heterosexual men", # Native American / Alaska Native
               "NAAN gay men",
               "NAAN bisexual men",
               "NAAN heterosexual women",
               "NAAN gay/lesbian women",
               "NAAN bisexual women",
               "NHPI heterosexual men", 
               "NHPI gay men", # Native Hawaiian / Pacific Islander
               "NHPI bisexual men",
               "NHPI heterosexual women",
               "NHPI gay/lesbian women",
               "NHPI bisexual women",
               "Multiracial heterosexual men",
               "Multiracial gay men",
               "Multiracial bisexual men",
               "Multiracial heterosexual women",
               "Multiracial gay/lesbian women",
               "Multiracial bisexual women")


# create a shorter version of the names
# used for intersectional group names in Figures 1 and 2
library(tm)
stopwords <- c(" women", " men")
groupNameShort <- tm::removeWords(groupName, 
                                  stopwords)
rm(stopwords) # delete this (not needed in further analysis)


# Autocorrelation plots
autocor_plots <- function(array, parsList) {
  lapply(parsList, function(x) {
    bayesplot::mcmc_acf(array, pars = x) +
      ggtitle(paste("Autocorrelation plots for", x)) +
      theme_clean() +
      theme(axis.title = element_text(face = "bold"),
            text = element_text(family = "Arial"))
  })
}


# Trace plots
trace_plots <- function(array, parsList) {
  lapply(parsList, function(x) {
    bayesplot::mcmc_trace(array, pars = x) +
      ggtitle(paste("Trace plots for", x)) +
      theme_clean() +
      theme(axis.title = element_text(face = "bold"),
            text = element_text(family = "Arial"))
  })
}


# Define logit2prob function
# This is used for converting logits to a 0-1 probability scale
# Input:
#    logit = a numeric value on the logit scale
logit2prob <- function(logit) {
  odds = exp(logit)
  prob = odds / (1 + odds)
  return(prob)
}



# Rename parameter values to match intersectional group names
# This is used to make the convergence diagnostics easier to read
# Input:
#    model = a stanfit model object
#    names = a vector of intersectional group names
#    start = start location (i.e., intercept for group 1)
#    finish = finish location (i.e., intercept for group 42)
renameArray <- function(model, names, start, finish) {
  
  # save model fit into an array object
  array <- as.array(model)
    
  # change names of array parameters
  for (i in start:finish) {
    dimnames(array)[[3]][i] <- names[i - start + 1]
  }
  
  # return the array with renamed parameters
  return(array)
}



# Define a function to calculate age-adjusted prevalence estimates
#    for each intersectional group. Age adjustment is calculated by
#    weighting the gamma parameter samples (fixed/main effects) for each age 
#    category to their population distribution. The weighted gammas are then 
#    combined with the random intercepts for each intersectional group and then 
#    the estimate is converted from logit to probability scale using the
#    `logit2prob` function defined above.
#
# Inputs:
#    model       = stanfit object
#    parameters  = parameter names from the stanfit object
#    age_weights = weights for each age category based on overall population
#                  distribution
#    return_samples = When FALSE (default), return summary statistics.
#                     When TRUE, return posterior samples.
age_adj_prevalence <- function(model, parameters, age_weights, 
                               return_samples = FALSE) {
  
  require(rstan)
  require(dplyr)
  
  ###############################################
  ### Extract main effects (gamma parameters) ###
  ###############################################
  
  # Main (fixed) effects: intercept & age category estimates
  # Intercept represents mean group value when age category = 18-25 (reference)
  gammas <- rstan::extract(model,
                           pars = parameters) %>% 
    as.data.frame() # convert to data frame
  
  # Rename parameters
  # Intercept is renamed "age1825"
  colnames(gammas)[1:4] <- 
    c("age1825", "age2634", "age3549", "age50plus")
  
  # Add the intercept (age1825) to the other age category estimates.
  # This puts each age category on the same relative scale, since the
  # age2634, age3549, and age50plus are only interpretable with the
  # reference value (age1825, the intercept).
  gammas <- gammas %>% 
    mutate(age2634 = age1825 + age2634,
           age3549 = age1825 + age3549,
           age50plus = age1825 + age50plus)
  
  
  ##################################################
  ### Extract random intercepts (u_j parameters) ###
  ##################################################
  
  # grab the parameter names for the random intercepts
  #    cs_sampling fits have "r_1_1[...]" format
  #    brms fits have "r_strata[...]" format
  pars_list <- grep("^r_1_1\\[\\d+\\]$|r_strata\\[\\d+,Intercept\\]", 
                    names(model), 
                    value = TRUE)
  
  # Extract random intercept samples for each intersectional group
  intercepts <- rstan::extract(model, pars = pars_list) %>% 
    as.data.frame() # convert to data frame
  
  # create an empty matrix to store prevalence estimates for each
  # intersectional group
  empty_matrix <- function(){ # define function to automate this process
    mat <- 
      matrix(data = NA,
             nrow = nrow(intercepts),
             ncol = ncol(intercepts))
    return(mat)
  } 
  age1825_samples <- empty_matrix()
  age2634_samples <- empty_matrix()
  age3549_samples <- empty_matrix()
  age50plus_samples <- empty_matrix()
  
  # loop through columns in the `intercepts` matrix
  for (i in 1:ncol(intercepts)) {
    # For each column (which represents an intersectional group),
    # sum the random intercept samples and age distribution-weighted gammas.
    # Then, transpose the samples so:
    #    rows = intersectional groups
    #    cols = posterior samples for each group
    age1825_samples[,i] <- intercepts[,i] + gammas$age1825
    age2634_samples[,i] <- intercepts[,i] + gammas$age2634 
    age3549_samples[,i] <- intercepts[,i] + gammas$age3549 
    age50plus_samples[,i] <- intercepts[,i] + gammas$age50plus 
  }
  
  # Transpose the matrices so that:
  #    rows = intersectional groups (42)
  #    cols = individual posterior samples (8000)
  # Then convert from logits to probability 0-1 scale
  age1825_samples <- age1825_samples %>% 
    t() %>% logit2prob()
  age2634_samples <- age2634_samples %>% 
    t() %>% logit2prob()
  age3549_samples <- age3549_samples %>% 
    t() %>% logit2prob()
  age50plus_samples <- age50plus_samples %>% 
    t() %>% logit2prob()
  
  # Create a unified matrix that adjusts for age category based on the
  # population-level age distribution. Weighted average of prevalence estimates
  # from each age category and intersectional group.
  adjusted_samples <- 
    age1825_samples*age_weights[1] + 
    age2634_samples*age_weights[2] +
    age3549_samples*age_weights[3] + 
    age50plus_samples*age_weights[4]
  
  # If needed, return the entire matrix of age-standardized posterior samples
  if (return_samples == TRUE) {
    return(adjusted_samples)
  }
  
  # Otherwise, create summary statistics of the posterior samples for
  # each intersectional group
  if (return_samples == FALSE) {
    # Create summary stats (mean + 95% CI) for each prevalence estiamte
    prevalence_summary <- 
      apply(adjusted_samples, 1,
            function(row) 
              # for each row (intersectional group), calculate mean and 95% CI
              c(mean(row),
                quantile(row, 
                         probs = c(0.025,0.975))))
    
    # transpose: rows = intersectional groups; cols = estimates
    # convert to data frame object
    prevalence_summary <- 
      t(prevalence_summary) %>% 
      as.data.frame() 
    
    # rounded to 1 digit on 0-100% scale
    prevalence_summary <- format(round(100*prevalence_summary,1), nsmall=1)
    
    # Rename the columns
    colnames(prevalence_summary) <- c("Mean", "LowerCI", "UpperCI")
    
    # Create a combined mean + 95% CI estimate column
    prevalence_summary <- prevalence_summary %>% 
      mutate(
        Combined = 
          paste(Mean, " (", LowerCI, ", ", UpperCI, ")", sep=""),
        # Remove extra white space
        Combined = gsub("\\(\\s+", "(", Combined)) %>% 
      # convert individual estimates back to numeric format
      mutate(across(c("Mean", "LowerCI", "UpperCI"), as.numeric)) %>% 
      select(Mean, LowerCI, UpperCI, Combined)
    
    
    # Final processing
    prevalence_summary <- prevalence_summary %>% 
      # add intersectional group names and IDs
      cbind(groupName, GroupID = c(1:nrow(.)), .) %>% 
      dplyr::rename("Intersectional Group" = groupName) %>% 
      mutate(
        # generate identity indicators (used for sorting)
        race = factor(
          rep(c(1:7),
              times=c(6,6,6,6,6,6,6)),
          labels = c("White", "Black", "Hispanic/Latine",
                     "Asian", "NAAN", "NHPI", "Multiracial")),
        gender = factor(
          rep(c(1,1,1,2,2,2),
              times=7), 
          labels = c("Men", "Women")),
        sexualOrientation = factor(
          rep(c(1,2,3),
              times=14),
          labels = c("Heterosexual", "Gay", "Bisexual"))) %>% 
      # Add shortened names for plotting Figures 1 & 2
      cbind(., groupNameShort)
    
    # output the prevalence estimates for each intersectional group
    return(prevalence_summary)
  }
  
}



# Define function to calculate SD, variance, VPC, and PCV
# VPC = variance partition coefficient
# PCV = proportional change in variance (relative to Model 1)
# Input:
#     models = a list of stanfit model objects
variance <- function(models) {
  
  # create a matrix to store mean + 95% CI estimates
  vpc <- matrix(nrow = length(models), ncol = 3)
  # define row/column names for clarity
  colnames(vpc) <- c("Mean VPC", "Lower CI", "Upper CI")
  rownames(vpc) < c("M1", "M2")
  
  # create matrix to store SD and variance values
  var_sd <- matrix(nrow = length(models), ncol = 2)
  # define row/column names for clarity
  colnames(var_sd) <- c("Variance", "SD")
  rownames(var_sd) <- c("M1", "M2")
  
  # loop through each model fit to extract/summarize var, SD, & VPC values
  for (i in 1:length(models)) {
    # using one stanfit object at a time
    model <- models[[i]]
    # extract posterior distribution of model standard deviation (SD)
    sd <- rstan::extract(model, pars = "sd_1")
    sd <- sd$sd_1
    # convert from SD to variance (sigma2)
    sigma2 <- sd^2
    
    # store variance and SD means
    var_sd[i,1] <- mean(sigma2)
    var_sd[i,2] <- mean(sd)
    
    # calculate full posterior distribution of VPC values
    fullPosterior <- sigma2/(sigma2+((pi^2)/3))*100
    # summarize VPC as mean and 95% CI (rounded to 2 decimal places)
    vpc[i,1] <- round(mean(fullPosterior),2)
    vpc[i,2:3] <- round(quantile(fullPosterior, 
                                prob = c(0.025, 0.975)),2)
  }
  
  # calculate proportional change in variance (PCV) relative to Model 1
  pcv <- matrix(nrow = nrow(vpc), ncol = 1)
  # define row/column names for clarity
  colnames(pcv) <- c("PCV")
  rownames(pcv) <- c("M1", "M2")
  # loop through Variance values to calculate PCV from Model 1 to 2
  for (j in 2:nrow(var_sd)) {
    pcv[j,1] <- round(((var_sd[1,1] - var_sd[j,1]) / var_sd[1,1])*100,2)
  }
  
  # save output into a list
  output <- list(var_sd, vpc, pcv)
  names(output) <- c("var_sd", "vpc", "pcv")
  
  # return the output file
  return(output)
}



# Define a function to calculate predicted prevalence from the full models
# based on total effects, main effects, and intersectional effects.
# Input:
#     model = a stanfit model object
#     parameters = a vector of main effects parameter names
#     age_weights = a vector of US adult population age distribution weights
calculate_effects <- function(model, parameters, age_weights) {
  
  
  ###############################################
  ### Extract main effects (gamma parameters) ###
  ###############################################
  
  # main/fixed effects: intercept + age category estimates
  # reference group (intercept) is age category 18-25
  gammas <- rstan::extract(model,
                           pars = parameters) %>% 
    as.data.frame() # convert to data frame
  
  # Rename variables
  # Call the intercept `age1825`
  colnames(gammas)[1:13] <- 
    c("age1825", "age2634", "age3549", "age50plus",
      "asian", "black", "latine", "naan", "nhpi", "multiracial",
      "woman", "gay", "bisexual")
  
  # add the intercept (age1825 reference) to the other age category estimates
  # this puts each age category on the same relative scale
  gammas <- gammas %>% 
    mutate(age2634 = age1825 + age2634,
           age3549 = age1825 + age3549,
           age50plus = age1825 + age50plus)
  
  
  #################################################
  ### Extract random intercepts (u_j parameters) ###
  ##################################################
  
  # grab the parameter names for the random intercepts
  #    cs_sampling fits have "r_1_1[...]" format
  #    brms fits have "r_strata[...]" format
  pars_list <- grep("^r_1_1\\[\\d+\\]$|r_strata\\[\\d+,Intercept\\]", 
                    names(model), 
                    value = TRUE)
  
  # Extract random intercept samples for each intersectional group
  intercepts <- rstan::extract(model, pars = pars_list) %>% 
    as.data.frame() %>% # convert to data frame
    t() # transpose so that rows = groups; cols = samples
  
  # create an empty matrix to store prevalence estimates for each
  # intersectional group
  empty_matrix <- function() { # define function to automate this process
    mat <- 
      matrix(data = NA,
             nrow = nrow(intercepts),
             ncol = ncol(intercepts))
    return(mat)
  } 
  
  age1825_samples <- empty_matrix()
  age2634_samples <- empty_matrix()
  age3549_samples <- empty_matrix()
  age50plus_samples <- empty_matrix()
  
  # loop through rows (i.e., each intersectional group)
  for (i in 1:nrow(intercepts)) {
    # For each row (which represents an intersectional group),
    # add the age distribution-weighted gammas.
    age1825_samples[i,] <- gammas$age1825
    age2634_samples[i,] <- gammas$age2634 
    age3549_samples[i,] <- gammas$age3549 
    age50plus_samples[i,] <- gammas$age50plus 
  }
  
  
  # define function to calculate main effects for each intersectional group
  # Input:
  #     samples = a matrix of posterior distribution samples where:
  #                   rows = intersectional groups (42 groups)
  #                   cols = posterior samples (8000 iterations)
  add_main_effects <- function(samples) {
    
    # Create an empty matrix to store main effect estimates
    main_effect_samples <- matrix(data = NA, nrow=42, ncol=8000)
    
    # Loop through each row
    for (i in 1:nrow(samples)) {
      
      # Define the first row as the baseline main effect
      # This represents White heterosexual men, where all predictors (except
      # age category) are = 0.
      main_effect_samples[i,] <- samples[i,]
      
      # Black main effect
      if (i %in% 7:12) {
        main_effect_samples[i,] = main_effect_samples[i,] + gammas$black
      }
      
      # Hispanic/Latine main effect
      if (i %in% 13:18) {
        main_effect_samples[i,] = main_effect_samples[i,] + gammas$latine
      }
      
      # Asian main effect
      if (i %in% 19:24) {
        main_effect_samples[i,] = main_effect_samples[i,] + gammas$asian
      }
      
      # Native American & Alaska Native main effect
      if (i %in% 25:30) {
        main_effect_samples[i,] = main_effect_samples[i,] + gammas$naan
      }
      
      # Native Hawaiian & Pacific Islander main effect
      if (i %in% 31:36) {
        main_effect_samples[i,] = main_effect_samples[i,] + gammas$nhpi
      }
      
      # Multiracial main effect
      if (i %in% 37:42) {
        main_effect_samples[i,] = main_effect_samples[i,] + gammas$multiracial
      }
      
      # Woman main effect
      if (i %in% c(4:6,10:12,16:18,22:24,28:30,34:36,40:42)) {
        main_effect_samples[i,] = main_effect_samples[i,] + gammas$woman
      }
      
      # Gay/lesbian main effect
      if (i %in% c(2,5,8,11,14,17,20,23,26,29,32,35,38,41)) {
        main_effect_samples[i,] = main_effect_samples[i,] + gammas$gay
      }
      
      # Bisexual main effect
      if (i %in% c(3,6,9,12,15,18,21,24,27,30,33,36,39,42)) {
        main_effect_samples[i,] = main_effect_samples[i,] + gammas$bisexual
      }
      
    }
    
    # return the main effect samples
    return(main_effect_samples)
  }
  
  # Add main effects for each age category
  # function is defined separately
  main_effects_1825 <- add_main_effects(age1825_samples)
  main_effects_2634 <- add_main_effects(age2634_samples)
  main_effects_3549 <- add_main_effects(age3549_samples)
  main_effects_50plus <- add_main_effects(age50plus_samples)
  
  
  # create empty matrices for total effects for each age category
  # total effects = main effects + random effects (intercepts)
  total_effects_1825 <- empty_matrix()
  total_effects_2634 <- empty_matrix()
  total_effects_3549 <- empty_matrix()
  total_effects_50plus <- empty_matrix()
  
  # loop through rows in the `intercepts` matrix
  for (i in 1:nrow(intercepts)) {
    # For each row (which represents an intersectional group),
    # sum the random intercept samples and the main effects
    total_effects_1825[i,] <- intercepts[i,] + main_effects_1825[i,]
    total_effects_2634[i,] <- intercepts[i,] + main_effects_2634[i,] 
    total_effects_3549[i,] <- intercepts[i,] + main_effects_3549[i,] 
    total_effects_50plus[i,] <- intercepts[i,] + main_effects_50plus[i,]
  }
  
  # summarize total effects
  total_effects <- empty_matrix() # create empty matrix to store values
  # for each age category, convert from logit to probability scale, then
  # multiply by the population age distribution weights, then sum element-wise, 
  # then multiply by 100 to get to a 0-100% scale
  total_effects <- 100 *
    (logit2prob(total_effects_1825) * age_weights[1] +
       logit2prob(total_effects_2634) * age_weights[2] +
       logit2prob(total_effects_3549) * age_weights[3] +
       logit2prob(total_effects_50plus) * age_weights[4])
  
  
  # summarize main effects
  main_effects <- empty_matrix() # create empty matrix to store values
  # for each age category, convert from logit to probability scale, then
  # multiply by the population age distribution weights, then sum element-wise, 
  # then multiply by 100 to get to a 0-100% scale
  main_effects <- 100 *
    (logit2prob(main_effects_1825) * age_weights[1] +
       logit2prob(main_effects_2634) * age_weights[2] +
       logit2prob(main_effects_3549) * age_weights[3] +
       logit2prob(main_effects_50plus) * age_weights[4]) 
  
  # summarize intersectional interaction effects
  intersectional_effects <- total_effects - main_effects # point difference
  intersectional_effects_pct <- # percentage difference
    100 * (total_effects - main_effects) / total_effects
  
  
  # create summaries
  # `summarize_effects` is a separately defined function
  total_sum <- summarize_effects(total_effects)
  main_sum <- summarize_effects(main_effects)
  intersectional_sum <- summarize_effects(intersectional_effects)
  intersectional_pct_sum <- summarize_effects(intersectional_effects_pct)
  
  
  # Return the following objects created in the function
  return(list(total_effects = total_effects, 
              main_effects = main_effects,
              intersectional_effects = intersectional_effects,
              total_sum = total_sum,
              main_sum = main_sum,
              intersectional_sum = intersectional_sum,
              intersectional_pct_sum = intersectional_pct_sum))
  
  
}


# Summarize the posterior distribution for the total effects, main effects, and
# intersectional (interaction) effects.
# Input:
#     input = posterior samples where rows are intersectional groups (42 groups)
#             and cols are posterior samples (8000 iterations)
summarize_effects <- function(input) {
  # Calculate mean, 2.5% quantile, and 97.5% quantile for each row
  # Each row is the set of posterior samples for an intersectional group
  summary_stats <- apply(input, 1, function(row) {
    # round values to 1 decimal place
    mean_val <- round(mean(row), 1)
    quantile_025 <- round(quantile(row, probs = 0.025), 1)
    quantile_975 <- round(quantile(row, probs = 0.975), 1)
    return(c(Mean = mean_val, 
             Quantile_025 = quantile_025, 
             Quantile_975 = quantile_975))
  })
  
  # Convert the result into a data frame
  summary_df <- as.data.frame(t(summary_stats)) 
  
  # Add intersectional group names to the rows
  rownames(summary_df) <- groupName
  
  # Edit column names 
  colnames(summary_df) <- c("Mean", "Lower CI", "Upper CI")
  
  # Create a "Combined" column that combines the mean and lower/upper CI
  # remove extra white space and ensure trailing zeros are preserved
  summary_df <- summary_df %>% 
    mutate(Combined = 
             paste(sprintf("%.1f", Mean), " (", 
                   gsub("(\\.\\d)0+(\\s*[,\\)])", "\\1\\2", 
                        sprintf("%.1f", `Lower CI`)), 
                   ", ", 
                   gsub("(\\.\\d)0+(\\s*[,\\)])", "\\1\\2", 
                        sprintf("%.1f", `Upper CI`)), 
                   ")", 
                   sep = "")
    )
  
  # Return the summary data frame
  return(summary_df)
}

