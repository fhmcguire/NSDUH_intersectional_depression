## ----------------------------------------------------------------------------
##
## Script name:       nsduh_analysis.R
##
## Purpose of script: Intersectional MAIHDA of past-year and lifetime major 
##                    depressive episode (MDE) among US adults at the 
##                    intersection of race/ethnicity, gender, and sexual 
##                    orientation.
##
## Author:            F. Hunter McGuire, MPH
##
## Date created:      2022-08-27
## Last update:       2023-04-06
##
## --------------------------------
##
## Data source: National Survey on Drug Use and Health (NSDUH) 2015-2020
##
## Input data: nsduh.rds
##
## Output data: Model fits (pyMDE_m1, pyMDE_m2a, pyMDE_m2b, pyMDE_m2c, pyMDE_m3)
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
##    'sp_1'       = subpop for past-year MDE analysis (0 = no, 1 = yes)
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
RMDpath = gsub('/02_nsduh_analysis.Rmd','', RMDpath)
### set the working directory
setwd(RMDpath)

# Load packages
library(csSampling)
library(parallel)

# Load data (cleaned in the "nsduh_dataManagement.R" file)
nsduh <- readRDS(file="nsduh.RDS")



#################
# PAST-YEAR MDE #
#################

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
# Create subpop-defined survey design object
# sp_2 = past-year MDE subpopulation
design_sp2 <- subset(design, sp_2==1) 
# Check mean of weights from subpop survey design object
summary(weights(design_sp2))
# Calibrate weights so the mean=1
design_sp2 <- calibrate(design_sp2,
                        ~ 1,
                        c(`(Intercept)` = length(weights(design_sp2))))
# Check mean after calibration
summary(weights(design_sp2)) # Mean=1


########################################
###   Define cs_sampling functions   ###
### Written by Matthew Williams, PhD ###
########################################

### Define new functions below then use above
list_2D_row_subset_DROP <- function (nmlist, rindex) 
{
  temp_list <- list()
  for (k in 1:length(nmlist)) {
    tmpdim <- dim(nmlist[[k]])
    ldim <- length(tmpdim)
    lcommas <- paste(rep(",", ldim - 1), collapse = " ")
    #copy over to new list - drop = FALSE retains ALL dimensions
    eval(parse(text = paste("temp_list$", 
                            names(nmlist)[k], 
                            " <- ", 
                            "(nmlist$", 
                            names(nmlist)[k], 
                            ")[rindex", 
                            lcommas, 
                            ",drop = FALSE]", 
                            sep = "")))
    #drop only first dimension of array - not others of size 1
    if(ldim > 1){
      eval(parse(text = paste("temp_list$", 
                              names(nmlist)[k], 
                              " <- ", 
                              "array(temp_list$", 
                              names(nmlist)[k], 
                              ", dim = tmpdim[-1])", 
                              sep = "")))
    }
    #if only had 1 dim which is the MCMC draw, make a double (no dim), 
    # rather than an array of dim 1 or 0
    if(ldim == 1){
      eval(parse(text = paste("temp_list$", 
                              names(nmlist)[k], 
                              " <- ", 
                              "as.double(temp_list$", 
                              names(nmlist)[k], 
                              ")", 
                              sep = "")))
    }	  	
  }
  return(temp_list)
}

### Define 'cs_sampling_DROP' function
cs_sampling_DROP <- function (svydes, 
                              mod_stan, 
                              par_stan = NA, 
                              data_stan, 
                              ctrl_stan = list(chains = 1, 
                                               iter = 2000, 
                                               warmup = 1000, 
                                               thin = 1), 
                              rep_design = FALSE, 
                              ctrl_rep = list(replicates = 100, 
                                              type = "mrbbootstrap"), 
                              sampling_args = list()) 
{
  require(rstan)
  require(survey)
  require(plyr)
  require(pkgcond)
  if (rep_design) {
    svyweights <- svydes$pweights
  }else {
    svyweights <- weights(svydes)
  }
  if (is.null(svyweights)) {
    if (!is.null(weights(data_stan))) {
      stop("No survey weights")
    }
  }
  if (is.null(weights(data_stan))) {
    if (!is.null(svyweights)) {
      warning("No stan data weights, using survey weights instead")
      data_stan$weights = weights(svydes)
    }
  }
  if (!isTRUE(all.equal(as.numeric(weights(data_stan)), 
                        as.numeric(svyweights)))) {
    stop("Survey weights and stan data weights do not match")
  }
  if (mean(weights(data_stan)) != 1) {
    stop("Mean of the weights is not 1")
  }
  print("stan fitting")
  out_stan <- do.call(sampling, c(list(object = mod_stan, 
                                       data = data_stan, 
                                       pars = par_stan, 
                                       chains = ctrl_stan$chains, 
                                       iter = ctrl_stan$iter, 
                                       warmup = ctrl_stan$warmup, 
                                       thin = ctrl_stan$thin), 
                                  sampling_args))
  par_samps_list <- rstan::extract(out_stan, 
                                   permuted = TRUE)
  if (anyNA(par_stan)) {
    par_stan <- names(par_samps_list)[-length(names(par_samps_list))]
  }
  par_samps <- as.matrix(out_stan, 
                         pars = par_stan)
  for (i in 1:dim(par_samps)[1]) {
    if (i == 1) {
      upar_samps <- unconstrain_pars(out_stan, 
                                     list_2D_row_subset_DROP(
                                       par_samps_list, i))
    }
    else {
      upar_samps <- rbind(upar_samps, 
                          unconstrain_pars(out_stan, 
                                           list_2D_row_subset_DROP(
                                             par_samps_list, i)))
    }
  }
  row.names(upar_samps) <- 1:dim(par_samps)[1]
  upar_hat <- colMeans(upar_samps)
  Hhat <- -1 * optimHess(upar_hat, 
                         gr = function(x) {
    grad_log_prob(out_stan, x)
  })
  if (rep_design == TRUE) {
    svyrep <- svydes
  }
  else {
    svyrep <- as.svrepdesign(design = svydes, 
                             type = ctrl_rep$type, 
                             replicates = ctrl_rep$replicates)
  }
  print("gradient evaluation")
  rep_tmp <- withReplicates(design = svyrep, 
                            theta = grad_par, 
                            stanmod = mod_stan, 
                            standata = data_stan, 
                            par_hat = upar_hat)
  Jhat <- vcov(rep_tmp)
  Hi <- solve(Hhat)
  V1 <- Hi %*% Jhat %*% Hi
  R1 <- chol(V1, pivot = TRUE)
  pivot <- attr(R1, "pivot")
  R1 <- R1[, order(pivot)]
  R2 <- chol(Hi, pivot = TRUE)
  pivot2 <- attr(R2, "pivot")
  R2 <- R2[, order(pivot2)]
  R2i <- solve(R2)
  R2iR1 <- R2i %*% R1
  upar_adj <- aaply(upar_samps, 
                    1, DEadj, 
                    par_hat = upar_hat, 
                    R2R1 = R2iR1, 
                    .drop = TRUE)
  for (i in 1:dim(upar_adj)[1]) {
    if (i == 1) {
      par_adj <- unlist(constrain_pars(out_stan, 
                                       upar_adj[i,])[par_stan])
    }
    else {
      par_adj <- rbind(par_adj, 
                       unlist(constrain_pars(out_stan,
                                             upar_adj[i, ])[par_stan]))
    }
  }
  row.names(par_adj) <- 1:dim(par_samps)[1]
  colnames(par_samps) <- colnames(par_adj)
  rtn = list(stan_fit = out_stan, 
             sampled_parms = par_samps, 
             adjusted_parms = par_adj)
  class(rtn) = c("cs_sampling", class(rtn))
  return(rtn)
}


cs_sampling_DROP_eigen <- function (svydes, 
                                    mod_stan, 
                                    par_stan = NA, 
                                    data_stan, 
                                    ctrl_stan = list(chains = 1,
                                                     ctrl_rep = list(
                                                       replicates = 100, 
                                                       type = "mrbbootstrap"), 
                                                     sampling_args = list())) 
{
  require(rstan)
  require(survey)
  require(plyr)
  require(pkgcond)
  if (rep_design) {
    svyweights <- svydes$pweights
  }else {
    svyweights <- weights(svydes)
  }
  if (is.null(svyweights)) {
    if (!is.null(weights(data_stan))) {
      stop("No survey weights")
    }
  }
  if (is.null(weights(data_stan))) {
    if (!is.null(svyweights)) {
      warning("No stan data weights, using survey weights instead")
      data_stan$weights = weights(svydes)
    }
  }
  if (!isTRUE(all.equal(as.numeric(weights(data_stan)), 
                        as.numeric(svyweights)))) {
    stop("Survey weights and stan data weights do not match")
  }
  if (mean(weights(data_stan)) != 1) {
    stop("Mean of the weights is not 1")
  }
  print("stan fitting")
  out_stan <- do.call(sampling, c(list(object = mod_stan, 
                                       data = data_stan, 
                                       pars = par_stan, 
                                       chains = ctrl_stan$chains, 
                                       iter = ctrl_stan$iter, 
                                       warmup = ctrl_stan$warmup, 
                                       thin = ctrl_stan$thin), 
                                  sampling_args))
  par_samps_list <- rstan::extract(out_stan, 
                                   permuted = TRUE)
  if (anyNA(par_stan)) {
    par_stan <- names(par_samps_list)[-length(names(par_samps_list))]
  }
  par_samps <- as.matrix(out_stan, 
                         pars = par_stan)
  for (i in 1:dim(par_samps)[1]) {
    if (i == 1) {
      upar_samps <- unconstrain_pars(out_stan, 
                                     list_2D_row_subset_DROP(
                                       par_samps_list, i))
    }
    else {
      upar_samps <- rbind(upar_samps, 
                          unconstrain_pars(out_stan,
                                           list_2D_row_subset_DROP(
                                             par_samps_list, i)))
    }
  }
  row.names(upar_samps) <- 1:dim(par_samps)[1]
  upar_hat <- colMeans(upar_samps)
  Hhat <- -1 * optimHess(upar_hat, gr = function(x) {
    grad_log_prob(out_stan, x)
  })
  if (rep_design == TRUE) {
    svyrep <- svydes
  }
  else {
    svyrep <- as.svrepdesign(design = svydes, 
                             type = ctrl_rep$type, 
                             replicates = ctrl_rep$replicates)
  }
  print("gradient evaluation")
  rep_tmp <- withReplicates(design = svyrep, 
                            theta = grad_par, 
                            stanmod = mod_stan, 
                            standata = data_stan, 
                            par_hat = upar_hat)
  Jhat <- vcov(rep_tmp)
  Hi <- solve(Hhat)
  V1 <- Hi %*% Jhat %*% Hi
  eigV <- eigen(V1, symmetric = TRUE)
  R1 <- diag(sqrt(abs(eigV$values)))%*%t(eigV$vectors) #small - zero values negative
  R2 <- chol(Hi)
  R2i <- solve(R2)
  R2i <- solve(R2)
  R2iR1 <- R2i %*% R1
  upar_adj <- aaply(upar_samps, 1, DEadj, par_hat = upar_hat, 
                    R2R1 = R2iR1, .drop = TRUE)
  for (i in 1:dim(upar_adj)[1]) {
    if (i == 1) {
      par_adj <- unlist(constrain_pars(out_stan, 
                                       upar_adj[i,])[par_stan])
    }
    else {
      par_adj <- rbind(par_adj, 
                       unlist(constrain_pars(out_stan,
                                             upar_adj[i, ])[par_stan]))
    }
  }
  row.names(par_adj) <- 1:dim(par_samps)[1]
  colnames(par_samps) <- colnames(par_adj)
  rtn = list(stan_fit = out_stan, sampled_parms = par_samps, 
             adjusted_parms = par_adj)
  class(rtn) = c("cs_sampling", class(rtn))
  return(rtn)
}



#####################################################
#### Model 1: Age-adjusted prevalence/null model ####
#####################################################

# Set up model and data
stancode_pyMDE_m1 <- make_stancode(brmsformula(
  as.factor(pyMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + (1 | strata)),
  data = nsduh,
  family = bernoulli(link="logit"),
  save_model = "brms_pyMDE_m1.stan")
modbrms_pyMDE_m1 <- stan_model("brms_pyMDE_m1.stan")
databrms_pyMDE_m1 <- make_standata(brmsformula(
  as.factor(pyMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + (1 | strata)),
  data = design_sp2$variables,
  family = bernoulli(link="logit"))

# Set design weights = Stan model weights
databrms_pyMDE_m1$weights <- weights(design_sp2)*length(weights(design_sp2))/sum(weights(design_sp2))

# Model estimation
set.seed(828) # for reproducible results
pyMDE_m1 <- cs_sampling_DROP(
  svydes = design_sp2, # use the subpop-defined survey object
  mod_stan = modbrms_pyMDE_m1,
  data_stan = databrms_pyMDE_m1,
  ctrl_stan = list(chains = 4,
                   iter = 5000, # 20,000 total iterations
                   warmup = 1000, # 4,000 warmups
                   thin = 1),
  sampling_args = list(cores = parallel::detectCores()))
# Save the model
saveRDS(pyMDE_m1, file="02_fits/pyMDE_m1.RDS")



#################################################
#### Model 2A: Race/ethnicity-adjusted model ####
#################################################

# Set up model and data
stancode_pyMDE_m2a <- make_stancode(brmsformula(
  as.factor(pyMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + 
    asian + black + hispanic + naan + nhpi + multi + (1 | strata)),
  data = nsduh,
  family = bernoulli(link="logit"),
  save_model = "brms_pyMDE_m2a.stan")
modbrms_pyMDE_m2a <- stan_model("brms_pyMDE_m2a.stan")
databrms_pyMDE_m2a <- make_standata(brmsformula(
  as.factor(pyMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + 
    asian + black + hispanic + naan + nhpi + multi + (1 | strata)),
  data = design_sp2$variables,
  family = bernoulli(link="logit"))

# Set design weights = Stan model weights
databrms_pyMDE_m2a$weights <- weights(design_sp2)*length(weights(design_sp2))/sum(weights(design_sp2))

# Model estimation
set.seed(828) # for reproducible results
pyMDE_m2a <- cs_sampling_DROP(
  svydes = design_sp2, # use the subpop-defined survey object
  mod_stan = modbrms_pyMDE_m2a,
  data_stan = databrms_pyMDE_m2a,
  ctrl_stan = list(chains = 4,
                   iter = 5000, # 20,000 total iterations
                   warmup = 1000, # 4,000 warmups
                   thin = 1),
  sampling_args = list(cores = parallel::detectCores()))
# Save the model
saveRDS(pyMDE_m2a, file="02_fits/pyMDE_m2a.RDS")



#########################################
#### Model 2B: Gender-adjusted model ####
#########################################

# Set up model and data
stancode_pyMDE_m2b <- make_stancode(brmsformula(
  as.factor(pyMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + woman + (1 | strata)),
  data = nsduh,
  family = bernoulli(link="logit"),
  save_model = "brms_pyMDE_m2b.stan")
modbrms_pyMDE_m2b <- stan_model("brms_pyMDE_m2b.stan")
databrms_pyMDE_m2b <- make_standata(brmsformula(
  as.factor(pyMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + woman + (1 | strata)),
  data = design_sp2$variables,
  family = bernoulli(link="logit"))

# Set design weights = Stan model weights
databrms_pyMDE_m2b$weights <- weights(design_sp2)*length(weights(design_sp2))/sum(weights(design_sp2))

# Model estimation
set.seed(828) # for reproducible results
pyMDE_m2b <- cs_sampling_DROP(
  svydes = design_sp2, # use the subpop-defined survey object
  mod_stan = modbrms_pyMDE_m2b,
  data_stan = databrms_pyMDE_m2b,
  ctrl_stan = list(chains = 4,
                   iter = 5000, # 20,000 total iterations
                   warmup = 1000, # 4,000 warmups
                   thin = 1),
  sampling_args = list(cores = parallel::detectCores()))
# Save the model
saveRDS(pyMDE_m2b, file="02_fits/pyMDE_m2b.RDS")



#####################################################
#### Model 2C: Sexual Orientation-adjusted model ####
#####################################################

# Set up model and data
stancode_pyMDE_m2c <- make_stancode(brmsformula(
  as.factor(pyMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + gay + bisexual + (1 | strata)),
  data = nsduh,
  family = bernoulli(link="logit"),
  save_model = "brms_pyMDE_m2c.stan")
modbrms_pyMDE_m2c <- stan_model("brms_pyMDE_m2c.stan")
databrms_pyMDE_m2c <- make_standata(brmsformula(
  as.factor(pyMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + gay + bisexual + (1 | strata)),
  data = design_sp2$variables,
  family = bernoulli(link="logit"))

# Set design weights = Stan model weights
databrms_pyMDE_m2c$weights <- weights(design_sp2)*length(weights(design_sp2))/sum(weights(design_sp2))

# Model estimation
set.seed(828) # for reproducible results
pyMDE_m2c <- cs_sampling_DROP(
  svydes = design_sp2, # use the subpop-defined survey object
  mod_stan = modbrms_pyMDE_m2c,
  data_stan = databrms_pyMDE_m2c,
  ctrl_stan = list(chains = 4,
                   iter = 5000, # 20,000 total iterations
                   warmup = 1000, # 4,000 warmups
                   thin = 1),
  sampling_args = list(cores = parallel::detectCores()))
# Save the model
saveRDS(pyMDE_m2c, file="02_fits/pyMDE_m2c.RDS")



##############################################################
#### Model 3: Intersectional interaction (fully adjusted) ####
##############################################################

# Set up model and data
stancode_pyMDE_m3 <- make_stancode(brmsformula(
  as.factor(pyMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus +
    asian + black + hispanic + naan + nhpi + multi + 
    woman + gay + bisexual + (1 | strata)),
  data = nsduh,
  family = bernoulli(link="logit"),
  save_model = "brms_pyMDE_m3.stan")
modbrms_pyMDE_m3 <- stan_model("brms_pyMDE_m3.stan")
databrms_pyMDE_m3 <- make_standata(brmsformula(
  as.factor(pyMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus +
    asian + black + hispanic + naan + nhpi + multi + 
    woman + gay + bisexual + (1 | strata)),
  data = design_sp2$variables,
  family = bernoulli(link="logit"))

# Set design weights = Stan model weights
databrms_pyMDE_m3$weights <- weights(design_sp2)*length(weights(design_sp2))/sum(weights(design_sp2))

# Model estimation
# Note: iterations/warmups and adapt_delta were increased 
#       to assist with model convergence issues.
set.seed(828) # for reproducible results
pyMDE_m3 <- cs_sampling_DROP(
  svydes = design_sp2, # use the subpop-defined survey object
  mod_stan = modbrms_pyMDE_m3,
  data_stan = databrms_pyMDE_m3,
  ctrl_stan = list(chains = 4,
                   iter = 10000, # 40,000 total iterations
                   warmup = 2500, # 10,000 warmups
                   thin = 1,
                   adapt_delta = 0.95), # increased from default (0.8)
  sampling_args = list(cores = parallel::detectCores()))
# Save the model
saveRDS(pyMDE_m3, file="02_fits/pyMDE_m3.RDS")





################
# LIFETIME MDE #
################

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
# Create subpop-defined survey design object
# sp_1 = lifetime MDE subpopulation
design_sp1 <- subset(design, sp_1==1) 
# Check mean of weights from subpop survey design object
summary(weights(design_sp1)) # Mean=1.18
# Calibrate weights so the mean=1
design_sp1 <- calibrate(design_sp1,
                        ~ 1,
                        c(`(Intercept)` = length(weights(design_sp1))))
# Check mean after calibration
summary(weights(design_sp1)) # Mean=1


################################################
#### Model 1: Age-adjusted prevalence model ####
################################################

# Set up model and data
stancode_ltMDE_m1 <- make_stancode(brmsformula(
  as.factor(ltMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + (1 | strata)),
  data = nsduh,
  family = bernoulli(link="logit"),
  save_model = "brms_ltMDE_m1.stan")
modbrms_ltMDE_m1 <- stan_model("brms_ltMDE_m1.stan")
databrms_ltMDE_m1 <- make_standata(brmsformula(
  as.factor(ltMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + (1 | strata)),
  data = design_sp1$variables,
  family = bernoulli(link="logit"))

# Set design weights = Stan model weights
databrms_ltMDE_m1$weights <- weights(design_sp1)*length(weights(design_sp1))/sum(weights(design_sp1))



# Model estimation
set.seed(828) # for reproducible results
ltMDE_m1 <- cs_sampling_DROP(
  svydes = design_sp1, # use the subpop-defined survey object
  mod_stan = modbrms_ltMDE_m1,
  data_stan = databrms_ltMDE_m1,
  ctrl_stan = list(chains = 4,
                   iter = 5000, # 20,000 total iterations
                   warmup = 1000, # 4,000 warmups
                   thin = 1),
  sampling_args = list(cores = parallel::detectCores()))
# Save the model
saveRDS(ltMDE_m1, file="02_fits/ltMDE_m1.RDS")



#################################################
#### Model 2A: Race/ethnicity-adjusted model ####
#################################################

# Set up model and data
stancode_ltMDE_m2a <- make_stancode(brmsformula(
  as.factor(ltMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + 
    asian + black + hispanic + naan + nhpi + multi + (1 | strata)),
  data = nsduh,
  family = bernoulli(link="logit"),
  save_model = "brms_ltMDE_m2a.stan")
modbrms_ltMDE_m2a <- stan_model("brms_ltMDE_m2a.stan")
databrms_ltMDE_m2a <- make_standata(brmsformula(
  as.factor(ltMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + 
    asian + black + hispanic + naan + nhpi + multi + (1 | strata)),
  data = design_sp1$variables,
  family = bernoulli(link="logit"))

# Set design weights = Stan model weights
databrms_ltMDE_m2a$weights <- weights(design_sp1)*length(weights(design_sp1))/sum(weights(design_sp1))

# Model estimation
set.seed(828) # for reproducible results
ltMDE_m2a <- cs_sampling_DROP(
  svydes = design_sp1, # use the subpop-defined survey object
  mod_stan = modbrms_ltMDE_m2a,
  data_stan = databrms_ltMDE_m2a,
  ctrl_stan = list(chains = 4,
                   iter = 5000, # 20,000 total iterations
                   warmup = 1000, # 4,000 warmups
                   thin = 1),
  sampling_args = list(cores = parallel::detectCores()))
# Save the model
saveRDS(ltMDE_m2a, file="02_fits/ltMDE_m2a.RDS")



#########################################
#### Model 2B: Gender-adjusted model ####
#########################################

# Set up model and data
stancode_ltMDE_m2b <- make_stancode(brmsformula(
  as.factor(ltMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + woman + (1 | strata)),
  data = nsduh,
  family = bernoulli(link="logit"),
  save_model = "brms_ltMDE_m2b.stan")
modbrms_ltMDE_m2b <- stan_model("brms_ltMDE_m2b.stan")
databrms_ltMDE_m2b <- make_standata(brmsformula(
  as.factor(ltMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + woman + (1 | strata)),
  data = design_sp1$variables,
  family = bernoulli(link="logit"))

# Set design weights = Stan model weights
databrms_ltMDE_m2b$weights <- weights(design_sp1)*length(weights(design_sp1))/sum(weights(design_sp1))

# Model estimation
set.seed(828) # for reproducible results
ltMDE_m2b <- cs_sampling_DROP(
  svydes = design_sp1, # use the subpop-defined survey object
  mod_stan = modbrms_ltMDE_m2b,
  data_stan = databrms_ltMDE_m2b,
  ctrl_stan = list(chains = 4,
                   iter = 5000, # 20,000 total iterations
                   warmup = 1000, # 4,000 warmups
                   thin = 1),
  sampling_args = list(cores = parallel::detectCores()))
# Save the model
saveRDS(ltMDE_m2b, file="02_fits/ltMDE_m2b.RDS")



#####################################################
#### Model 2C: Sexual Orientation-adjusted model ####
#####################################################

# Set up model and data
stancode_ltMDE_m2c <- make_stancode(brmsformula(
  as.factor(ltMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + gay + bisexual + (1 | strata)),
  data = nsduh,
  family = bernoulli(link="logit"),
  save_model = "brms_ltMDE_m2c.stan")
modbrms_ltMDE_m2c <- stan_model("brms_ltMDE_m2c.stan")
databrms_ltMDE_m2c <- make_standata(brmsformula(
  as.factor(ltMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus + gay + bisexual + (1 | strata)),
  data = design_sp1$variables,
  family = bernoulli(link="logit"))

# Set design weights = Stan model weights
databrms_ltMDE_m2c$weights <- weights(design_sp1)*length(weights(design_sp1))/sum(weights(design_sp1))

# Model estimation
set.seed(828) # for reproducible results
ltMDE_m2c <- cs_sampling_DROP(
  svydes = design_sp1, # use the subpop-defined survey object
  mod_stan = modbrms_ltMDE_m2c,
  data_stan = databrms_ltMDE_m2c,
  ctrl_stan = list(chains = 4,
                   iter = 5000, # 20,000 total iterations
                   warmup = 1000, # 4,000 warmups
                   thin = 1),
  sampling_args = list(cores = parallel::detectCores()))
# Save the model
saveRDS(ltMDE_m2c, file="02_fits/ltMDE_m2c.RDS")



##############################################################
#### Model 3: Intersectional interaction (fully adjusted) ####
##############################################################

# Set up model and data
stancode_ltMDE_m3 <- make_stancode(brmsformula(
  as.factor(ltMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus +
    asian + black + hispanic + naan + nhpi + multi + 
    woman + gay + bisexual + (1 | strata)),
  data = nsduh,
  family = bernoulli(link="logit"),
  save_model = "brms_ltMDE_m3.stan")
modbrms_ltMDE_m3 <- stan_model("brms_ltMDE_m3.stan")
databrms_ltMDE_m3 <- make_standata(brmsformula(
  as.factor(ltMDE)|weights(wt_norm) ~ 
    1 + age2634 + age3549 + age50plus +
    asian + black + hispanic + naan + nhpi + multi + 
    woman + gay + bisexual + (1 | strata)),
  data = design_sp1$variables,
  family = bernoulli(link="logit"))

# Set design weights = Stan model weights
databrms_ltMDE_m3$weights <- weights(design_sp1)*length(weights(design_sp1))/sum(weights(design_sp1))

# Model estimation
# Note: iterations/warmups and adapt_delta were increased to assist with model convergence issues
set.seed(828) # for reproducible results
ltMDE_m3 <- cs_sampling_DROP(
  svydes = design_sp1, # use the subpop-defined survey object
  mod_stan = modbrms_ltMDE_m3,
  data_stan = databrms_ltMDE_m3,
  ctrl_stan = list(chains = 4,
                   iter = 10000, # 40,000 total iterations
                   warmup = 2500, # 10,000 warmups
                   thin = 1,
                   adapt_delta = 0.95), # increased from default (0.8)
  sampling_args = list(cores = parallel::detectCores()))
# Save the model
saveRDS(ltMDE_m3, file="02_fits/ltMDE_m3.RDS")

