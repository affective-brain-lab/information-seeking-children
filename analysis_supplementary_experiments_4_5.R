# AUTHOR: GAIA MOLINARO
# Data analysis for supplementary experiments 4-5
# in Molinaro, Cogliati Dezza, Buehler, Moutsiana & Sharot (in prep.)

# CLEAR WORKSPACE
rm(list = ls())

# SET WORKING DIRECTORY
wd <- "insert_your_wd_here"
setwd(wd)

packages <- c("ggplot2", "GGally", "ggpubr", "psych", "nls2", "RColorBrewer",
              "grid", "optimx", "ggcorrplot", "tidyverse", "glue", "gridExtra", 
              "lme4", "table1", "BayesFactor", "interactions", "modelbased", "quest",
              "glmnet", "glmnetUtils", "broom", "MuMIn")
pacman::p_load(packages, character.only = TRUE)

# Additional basic functions
# Function to get the standard error from the mean
sem <- function(x) {
  return (sdamr::sample_sd(x)/sqrt(length(x)))
}

# COLOR SETTINGS
ev_color <- "#E41A1C"
unc_color <- "#377EB8"
agn_color <- "#4DAF4A"

# FUNCTIONS

# convert from raw to z
find_Z <- function(values, var) {
  return ((values - mean(var))/sd(var))
}

# convert from z to raw
find_non_Z <- function(values, var) {
  return ((values*sd(var)) + mean(var))
}


#### SUPPLEMENTARY EXPERIMENTS 4-5 (ADULTS) ####
#### Load and prep datasets #### 
rescale_variables <- TRUE  # whether to rescale the numeric variables
center_variables <- TRUE # whether to center numeric variables (except deltas)
rescale_age <- TRUE  # whether to rescale age
center_age <- TRUE # whether to center age

###### Load adults' data ###### 
data_batch <- "" # either 4 or 5
adu_dat <- read.csv(sprintf("data_supplementary_experiment_%s.csv", data_batch))
adu_dat <- subset(adu_dat, !(condition %in% c("catch_1", "catch_2")))
adu_dat <- subset(adu_dat, catch_trials_score == 100)

# Turn variables into appropriate data types

adu_dat$delta_EV_non_Z <- as.numeric(adu_dat$delta_EV)
adu_dat$delta_uncertainty_level_non_Z <- as.numeric(adu_dat$delta_uncertainty_level)
adu_dat$delta_agency_non_Z <- as.numeric(adu_dat$delta_agency)
adu_dat$delta_SD <- as.numeric(adu_dat$delta_SD)
adu_dat$age_in_years_non_Z <- as.numeric(adu_dat$age_in_years)
adu_dat$age_in_years <- scale(as.numeric(adu_dat$age_in_years), center=center_age, scale=rescale_age)
adu_dat$RT_info_choice_non_Z <- as.numeric(adu_dat$RT_info_choice)
adu_dat$RT_info_choice_log <- scale(log(as.numeric(adu_dat$RT_info_choice_non_Z)), center=center_variables, scale=rescale_variables)
adu_dat$RT_info_choice <- scale(as.numeric(adu_dat$RT_info_choice), center=center_variables, scale=rescale_variables)
adu_dat$subject_ID <- as.factor(adu_dat$subject_ID)
adu_dat$gender_coded[adu_dat$gender_coded==3] <- 1 # contrast will be female vs non-female
adu_dat$gender_coded <- as.factor(adu_dat$gender_coded)
contrasts(adu_dat$gender_coded) <- contr.helmert(2)
adu_dat$info_choice <- as.factor(adu_dat$info_choice) # 0 = left, 1 = right
adu_dat$percent_comprehension_non_Z <- as.numeric(adu_dat$percent_comprehension)
adu_dat$percent_comprehension <- scale(as.numeric(adu_dat$percent_comprehension), center=center_variables, scale=rescale_variables)
adu_dat$wob_non_Z <- as.numeric(adu_dat$wob)
adu_dat$wob <- scale(as.numeric(adu_dat$wob), center=center_variables, scale=rescale_variables)
adu_dat$fishing_rewardL <- as.numeric(adu_dat$fishing_rewardL)
adu_dat$fishing_rewardR <- as.numeric(adu_dat$fishing_rewardR)# correct fishing choices based on what participants see
adu_dat$fishing_non_Z <- as.numeric(adu_dat$prop_correct_fishing)
adu_dat$fishing <- scale(as.numeric(adu_dat$prop_correct_fishing), center=center_variables, scale=rescale_variables)

## Compute scaled deltas
# Scale without centering so the sign is preserved
rescale_deltas <- function(dat, deltas=c("EV", "uncertainty_level", "agency")) {
  if ("EV" %in% deltas){
    dat$delta_EV <- scale(dat$delta_EV_non_Z, center=FALSE)
  }
  if ("uncertainty_level" %in% deltas){
    dat$delta_uncertainty_level <- scale(dat$delta_uncertainty_level_non_Z, center=FALSE)
  }
  if ("agency" %in% deltas){
    dat$delta_agency <- scale(dat$delta_agency_non_Z, center=FALSE)
  }
  return(dat)
}

# Find min and max and rescale between -1 and 1
if (rescale_variables) {adu_dat <- rescale_deltas(adu_dat)}

#### Demographics ####
summary(subset(aud_dat, trial_number==1)$age_in_years_non_Z)
sdamr::sample_sd(subset(aud_dat, trial_number==1)$age_in_years_non_Z)
table(subset(dat, trial_number==1)$gender)

#### Comprehension, wob, and fishing scores ####
## t-tests
t.test(subset(adu_dat, trial_number == 1)$wob_non_Z, mu=0.5, var.equal=TRUE)

## Fishing choices
# t-tests

# Overall fishing choices compared to chance and children vs adults
# based on what participants see
t.test(subset(adu_dat, trial_number == 1)$fishing_non_Z, 
       mu=0.5, var.equal=TRUE)


#### Correlations between deltas ####
##### Matrices #####
adu_dat_corr <- adu_dat[corr_variables]
adu_dat_corr_mat <- round(cor(adu_dat_corr, method="pearson"), 3)
adu_dat_p_mat <- ggcorrplot::cor_pmat(adu_dat_corr)
ggcorrplot(adu_dat_corr_mat, outline.col = "white", lab = TRUE, insig = "blank")
adu_dat_corr_mat 
adu_dat_p_mat


####  Models of adults only #### 
{
  # Adults

  # Full 
  adu_mod_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + fishing + 
                                (delta_EV + delta_uncertainty_level + delta_agency 
                                 | subject_ID), 
                              data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without intercept
  adu_mod_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                     percent_comprehension + wob + fishing -1 +
                                                     (delta_EV + delta_uncertainty_level + delta_agency 
                                                      | subject_ID), 
                                                   data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without EV as a fixed effect
  adu_mod_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob + fishing + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  adu_mod_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                             percent_comprehension + wob + fishing + 
                                             (delta_EV + delta_uncertainty_level + delta_agency 
                                              | subject_ID), 
                                           data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  adu_mod_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                              percent_comprehension + wob + fishing + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without comprehension as a fixed effect
  adu_mod_full_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              wob + fishing +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without wob as a fixed effect
  adu_mod_full_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without fishing as a fixed effect
  adu_mod_full_drop_fishing_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                  percent_comprehension + wob  +
                                                  (delta_EV + delta_uncertainty_level + delta_agency
                                                   | subject_ID), 
                                                data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop EV random
  adu_mod_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob  + fishing +
                                               (delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  adu_mod_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob  + fishing +
                                              (delta_EV + delta_agency
                                               | subject_ID), 
                                            data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  adu_mod_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob  + fishing +
                                               (delta_EV + delta_uncertainty_level
                                                | subject_ID), 
                                             data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Affect only
  adu_mod_Af <- glmer(info_choice ~ delta_EV +
                              percent_comprehension + wob + fishing + 
                              (delta_EV 
                               | subject_ID), 
                            data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only
  adu_mod_C <- glmer(info_choice ~  delta_uncertainty_level +
                             percent_comprehension + wob + fishing + 
                             (delta_uncertainty_level 
                              | subject_ID), 
                           data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only 
  adu_mod_Ac <- glmer(info_choice ~ delta_agency +
                              percent_comprehension + wob + fishing +
                              (delta_agency 
                               | subject_ID), 
                            data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition 
  adu_mod_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               percent_comprehension + wob + fishing + 
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID), 
                             data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  adu_mod_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                percent_comprehension + wob + fishing + 
                                (delta_EV + delta_agency
                                 | subject_ID), 
                              data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  adu_mod_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + fishing + 
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID), 
                             data = adu_dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
}   

###### Significance tests ###### 
# access results of test a with a$Chisq, a$Df, a$Pr[2]
adu_test_fixed_EV <- anova(adu_mod_full, adu_mod_full_drop_Af_fixed, test="Chisq")
adu_test_fixed_uncertainty <- anova(adu_mod_full, adu_mod_full_drop_C_fixed, test="Chisq")
adu_test_fixed_agency <- anova(adu_mod_full, adu_mod_full_drop_Ac_fixed, test="Chisq")
# covariates
adu_test_fixed_pc <- anova(adu_mod_full, adu_mod_full_drop_pc_fixed, test="Chisq")
adu_test_fixed_wob <- anova(adu_mod_full, adu_mod_full_drop_wob_fixed, test="Chisq")
adu_test_fixed_fishing <- anova(adu_mod_full, adu_mod_full_drop_fishing_fixed, test="Chisq")
