# -----------------------------------------------------------------------------
#
# Script name: 01_nsduh_dataManagement.R
#
# Purpose of script: Data management for intersectional MAIHDA of lifetime
#                    major depressive episode (MDE) and past-year MDE among 
#                    US adults at the intersection of race/ethnicity, 
#                    sex/gender, and sexual orientation.
#
# Author: F. Hunter McGuire, PhD, MPH
#
# Date created: 2022-08-27
# Last update:  2023-10-11
#
# -----------------------------------------------------------------------------
#
# Input data:  NSDUH_2015.dta, NSDUH_2016.dta, NSDUH_2017.dta
#              NSDUH_2018.dta, NSDUH_2019.dta, NSDUH_2020.dta
#
# Output data: nsduh.rds (R format)
#
# Data source: National Survey on Drug Use and Health (NSDUH) 2015-2020 
#
# Variables:
#
#    'ltMDE'      = lifetime major depressive episode (0 = no, 1 = yes)
#    'pyMDE'      = past-year major depressive episode (0 = no, 1 = yes)
#    'mhservice'  = past-year mental health service use (0 = no, 1 = yes)
#
#    'asian'      = race indicator (1 = Asian)
#    'black'      = race indicator (1 = Black / African-American)
#    'hispanic'   = race indicator (1 = Hispanic)
#    'naan'       = race indicator (1 = Native American / Alaskan Native)
#    'nhpi'       = race indicator (1 = Native Hawaiian / Pacific Islander)
#    'multi'      = race indicator (1 = Multiracial)
#    'white'      = race indicator (1 = White)
#
#    'woman'      = gender indicator (1 = Woman)
#    'man'        = gender indicator (1 = Man)
#
#    'hetero'     = sexual orientation indicator (1 = Heterosexual)
#    'gay'        = sexual orientation indicator (1 = Gay / Lesbian)
#    'bisexual'   = sexual orientation indicator (1 = Bisexual)
#
#    'age1217'    = aged 12-17 years old (0 = no, 1 = yes)
#    'age1825'    = aged 18-25 years old (0 = no, 1 = yes)
#    'age2634'    = aged 26-34 years old (0 = no, 1 = yes)
#    'age3549'    = aged 35-49 years old (0 = no, 1 = yes)
#    'age50plus'  = aged 50+ years old (0 = no, 1 = yes)
#    
#    'sp_1'       = subpop for lifetime MDE analysis (0 = no, 1 = yes)
#    'sp_2'       = subpop for past-year MDE analysis (0 = no, 1 = yes)
#
#    'verep'      = NSDUH cluster weights
#    'vestr'      = NSDUH strata weights
#    'adjwt6'     = ADJUSTED NSDUH person-level weights for number of data 
#                   collection years (6 years: 2015-2020)
#
# -----------------------------------------------------------------------------



#################
##### SETUP #####
#################

# load packages
library(haven) # load Stata-format datasets
library(labelled) # manipulating labelled data
library(tidyverse) # general data management tasks
library(naniar) # recoding missing values

## set working directory to source file location
### get full file path
Rpath <- rstudioapi::getSourceEditorContext()$path
### remove file name from file path
Rpath <- gsub('/01_nsduh_dataManagement.R','', Rpath)
### set the working directory
setwd(Rpath)



#####################
##### LOAD DATA #####
#####################

# load full data files (all variables) for each NSDUH data collection year
nsduh2015 <- read_dta("01_original_data/NSDUH_2015.dta")
nsduh2016 <- read_dta("01_original_data/NSDUH_2016.dta")
nsduh2017 <- read_dta("01_original_data/NSDUH_2017.dta")
nsduh2018 <- read_dta("01_original_data/NSDUH_2018.dta")
nsduh2019 <- read_dta("01_original_data/NSDUH_2019.dta")
nsduh2020 <- read_dta("01_original_data/NSDUH_2020.dta")

# define function to extract necessary variables
extractVariables <- function(data) {
  
  # subset variables
  dataSubset <- data %>% 
    select(amdeyr, amdelt, amhtxrc3, 
           newrace2, irsex, sexident, catag6, 
           analwt_c, vestr, verep, year)
  # remove Stata formats and labels
  dataSubset <- remove_labels(dataSubset)
  dataSubset <- zap_formats(dataSubset)
  # return the smaller data file
  return(dataSubset)
  
}
# extract variables from each data set
nsduh2015new <- extractVariables(nsduh2015)
nsduh2016new <- extractVariables(nsduh2016)
nsduh2017new <- extractVariables(nsduh2017)
nsduh2018new <- extractVariables(nsduh2018)
nsduh2019new <- extractVariables(nsduh2019)
nsduh2020new <- extractVariables(nsduh2020)

# row bind data across NSDUH data collection years
nsduh <- rbind(nsduh2015new, nsduh2016new,
               nsduh2017new, nsduh2018new,
               nsduh2019new, nsduh2020new)



####################
#### MANAGEMENT ####
####################

# Recode values (refused, don't know, non-response) to missing 
nsduh <- nsduh %>% 
  naniar::replace_with_na(
    replace = list(
      sexident = c(85, 89, 94, 97, 98, 99)))

# General data management
nsduh <- nsduh %>% 
  # indicator variables for model variables
  mutate( 
    white = ifelse(newrace2 == 1, 1, 0), # Non-Hispanic/Latine (NHL) White
    black = ifelse(newrace2 == 2, 1, 0), # NHL Black
    naan = ifelse(newrace2 == 3, 1, 0), # NHL Native American
    nhpi = ifelse(newrace2 == 4, 1, 0), # NHL Native Hawaiian / Pac. Islander
    asian = ifelse(newrace2 == 5, 1, 0), # NHL Asian
    multi = ifelse(newrace2 == 6, 1, 0), # NHL Multiracial
    hispanic = ifelse(newrace2 == 7, 1, 0), # Hispanic/Latine
    man = ifelse(irsex == 1, 1, 0), # Man
    woman = ifelse(irsex == 2, 1, 0), # Woman
    hetero = ifelse(sexident == 1, 1, 0), # Heterosexual
    gay = ifelse(sexident == 2, 1, 0), # Gay / Lesbian
    bisexual = ifelse(sexident == 3, 1, 0), # Bisexual
    ltMDE = ifelse(amdelt == 1, 1, 0), # Lifetime MDE
    pyMDE = ifelse(amdeyr == 1, 1, 0), # Past-year MDE
    mhservice = ifelse(amhtxrc3 == 1, 1, 0), # Past-year MH service use
    age1217 = ifelse(catag6 == 1, 1, 0), # aged 12-17
    age1825 = ifelse(catag6 == 2, 1, 0), # aged 18-25
    age2634 = ifelse(catag6 == 3, 1, 0), # aged 26-34
    age3549 = ifelse(catag6 == 4, 1, 0), # aged 35-49
    age50plus = ifelse(catag6 == 5 | catag6 == 6, 1, 0) # aged 50+
    ) %>% 
  # subpopulation indicators
  mutate( 
    # lifetime MDE analysis subpopulation
    sp_1 = ifelse(!is.na(newrace2) & # non-missing data on model variables
                    !is.na(irsex) & 
                    !is.na(sexident) & 
                    !is.na(catag6) &
                    catag6 >= 2 & # 18 years or older
                    !is.na(amdelt), 1, 0),
    # past-year MDE analysis subpopulation
    sp_2 = ifelse(!is.na(newrace2) & # non-missing data on model variables
                    !is.na(irsex) & 
                    !is.na(sexident) & 
                    !is.na(catag6) &
                    catag6 >= 2 & # 18 years or older
                    !is.na(amdeyr), 1, 0))


# Create intersectional group variable (`strata`)
### Defined as mutually exclusive groups by race/ethnicity,
### gender, and sexual orientation.
nsduh <- nsduh %>%   
  mutate(
    strata = case_when(
      white == 1 & man == 1 & hetero == 1 ~ 1, # White heterosexual men
      white == 1 & man == 1 & gay == 1 ~ 2, # White gay men
      white == 1 & man == 1 & bisexual == 1 ~ 3, # White bisexual men
      white == 1 & woman == 1 & hetero == 1 ~ 4, # White heterosexual women
      white == 1 & woman == 1 & gay == 1 ~ 5, # White gay/lesbian women
      white == 1 & woman == 1 & bisexual == 1 ~ 6, # White bisexual women
      black == 1 & man == 1 & hetero == 1 ~ 7, # Black heterosexual men
      black == 1 & man == 1 & gay == 1 ~ 8, # Black gay men
      black == 1 & man == 1 & bisexual == 1 ~ 9, # Black bisexual men
      black == 1 & woman == 1 & hetero == 1 ~ 10, # Black heterosexual women
      black == 1 & woman == 1 & gay == 1 ~ 11, # Black gay/lesbian women
      black == 1 & woman == 1 & bisexual == 1 ~ 12, # Black bisexual women
      hispanic == 1 & man == 1 & hetero == 1 ~ 13, # Hispanic heterosexual men
      hispanic == 1 & man == 1 & gay == 1 ~ 14, # Hispanic gay men
      hispanic == 1 & man == 1 & bisexual == 1 ~ 15, # Hispanic bisexual men
      hispanic == 1 & woman == 1 & hetero == 1 ~ 16, # Hispanic heterosexual women
      hispanic == 1 & woman == 1 & gay == 1 ~ 17, # Hispanic gay/lesbian women
      hispanic == 1 & woman == 1 & bisexual == 1 ~ 18, # Hispanic bisexual women
      asian == 1 & man == 1 & hetero == 1 ~ 19, # Asian heterosexual men
      asian == 1 & man == 1 & gay == 1 ~ 20, # Asian gay men
      asian == 1 & man == 1 & bisexual == 1 ~ 21, # Asian bisexual men
      asian == 1 & woman == 1 & hetero == 1 ~ 22, # Asian heterosexual women
      asian == 1 & woman == 1 & gay == 1 ~ 23, # Asian gay/lesbian women
      asian == 1 & woman == 1 & bisexual == 1 ~ 24, # Asian bisexual women
      naan == 1 & man == 1 & hetero == 1 ~ 25, # Native American heterosexual men
      naan == 1 & man == 1 & gay == 1 ~ 26, # Native American gay men
      naan == 1 & man == 1 & bisexual == 1 ~ 27, # Native American bisexual men
      naan == 1 & woman == 1 & hetero == 1 ~ 28, # Native American heterosexual women
      naan == 1 & woman == 1 & gay == 1 ~ 29, # Native American gay/lesbian women
      naan == 1 & woman == 1 & bisexual == 1 ~ 30, # Native American bisexual women
      nhpi == 1 & man == 1 & hetero == 1 ~ 31, # Native Hawaiian / PI heterosexual men
      nhpi == 1 & man == 1 & gay == 1 ~ 32, # Native Hawaiian / PI gay men
      nhpi ==1 & man == 1 & bisexual == 1 ~ 33, # Native Hawaiian / PI bisexual men
      nhpi == 1 & woman ==1 & hetero == 1 ~ 34, # Native Hawaiian / PI heterosexual women
      nhpi == 1 & woman == 1 & gay == 1 ~ 35, # Native Hawaiian / PI gay/lesbian women
      nhpi == 1 & woman == 1 & bisexual == 1 ~ 36, # Native Hawaiian / PI bisexual women
      multi == 1 & man == 1 & hetero == 1 ~ 37, # Multiracial heterosexual men
      multi == 1 & man == 1 & gay == 1 ~ 38, # Multiracial gay men
      multi == 1 & man == 1 & bisexual == 1 ~ 39, # Multiracial bisexual men
      multi == 1 & woman == 1 & hetero == 1 ~ 40, # Multiracial heterosexual women
      multi == 1 & woman == 1 & gay == 1 ~ 41, # Multiracial gay/lesbian women
      multi == 1 & woman == 1 & bisexual == 1 ~ 42)) # Multiracial bisexual women

# Adjust the person-level NSDUH sampling weights for multiple years
# divide by the number of NSDUH data collection years (6 years: 2015-2020)
nsduh$adjwt6 <- nsduh$analwt_c/6 

# Select necessary variables for analysis
nsduh <- nsduh %>% 
  select(strata, # intersectional groups/strata
         asian, black, hispanic, naan, nhpi, multi, white, # race/ethnicity
         woman, man, # gender
         gay, bisexual, hetero, # sexual orientation
         age1825, age2634, age3549, age50plus, # age category
         ltMDE, pyMDE, mhservice, # outcome measures
         vestr, verep, adjwt6, # survey design weights
         year, # data collection year
         sp_1, sp_2) # subpopulation indicators

# Save the final data frame to be used in analysis
saveRDS(nsduh, file="nsduh.RDS")




# For Methods section on analytic sample:
# Calculate analytic sample exclusion criteria
samples <- nsduh %>% 
  # step 1: older than 18
  mutate(step1 = ifelse(age1825==0 & age2634==0 &
                          age3549==0 & age50plus==0, 0, 1),
         # step 2: non-missing sexual orientation
         step2 = ifelse(step1 == 1 & !is.na(hetero), 1, 0),
         # step 3: non-missing past-year MDE
         step3 = ifelse(step2 == 1 & !is.na(ltMDE), 1, 0),
         # step 4: non-missing lifetime MDE
         step4 = ifelse(step3 == 1 & !is.na(pyMDE), 1, 0))


table(samples$step1)
table(samples$step2)
table(samples$step3)

# missing sexual orientation
table(samples$step1)[2] - table(samples$step2)[2] 
# missing lifetime MDE
table(samples$step2)[2] - table(samples$step3)[2] 
# missing past-year MDE
table(samples$step3)[2] - table(samples$step4)[2] 
