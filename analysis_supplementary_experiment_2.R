# AUTHOR: GAIA MOLINARO
# Data analysis for supplementary experiment 2
# in Molinaro, Cogliati Dezza, & Sharot (in prep.)

# CLEAR WORKSPACE
rm(list = ls())

# SET WORKING DIRECTORY
wd = "insert_your_wd_here"

#### LOAD REQUIRED PACKAGES #####
packages <- c("psych", "optimx")
pacman::p_load(packages, character.only = TRUE)


#### SUPPLEMENTARY EXPERIMENT 2 (COGNITIVE UTILITY CONTROL) ####
#### SE2 Load and prep datasets #### 
setwd(wd)
dat_cog_ctrl <- read.csv("data_supplementary_experiment_2.csv")
dat_cog_ctrl <- subset(dat_cog_ctrl, !(condition %in% c("catch_1", "catch_2")))
dat_cog_ctrl <- subset(dat_cog_ctrl, !(age_in_years == 13))
dat_cog_ctrl <- subset(dat_cog_ctrl, catch_trials_score == 100)

# Create age groups
dat_cog_ctrl$age_group <- "None"
dat_cog_ctrl$age_group[dat_cog_ctrl$age_in_years  %in% c(4, 5)] <- "4-5"
dat_cog_ctrl$age_group[dat_cog_ctrl$age_in_years  %in% c(6, 7)] <- "6-7"
dat_cog_ctrl$age_group[dat_cog_ctrl$age_in_years  %in% c(8, 9)] <- "8-9"
dat_cog_ctrl$age_group[dat_cog_ctrl$age_in_years  %in% c(10, 11, 12)] <- "10-12"
dat_cog_ctrl$age_group <- factor(dat_cog_ctrl$age_group, levels=c("4-5", "6-7", "8-9","10-12"))
dat_cog_ctrl$age_group_coded <- "None"
dat_cog_ctrl$age_group_coded[dat_cog_ctrl$age_in_years  %in% c(4, 5)] <- 0
dat_cog_ctrl$age_group_coded[dat_cog_ctrl$age_in_years  %in% c(6, 7)] <- 1
dat_cog_ctrl$age_group_coded[dat_cog_ctrl$age_in_years  %in% c(8, 9)] <- 2
dat_cog_ctrl$age_group_coded[dat_cog_ctrl$age_in_years  %in% c(10, 11, 12)] <- 3
dat_cog_ctrl$age_group_coded <- factor(dat_cog_ctrl$age_group_coded, levels=c("4-5", "6-7", "8-9","10-12"))


# Turn variables into appropriate data types
# Use scale() for standardized values
dat_cog_ctrl$delta_uncertainty_level_non_Z <- (as.numeric(dat_cog_ctrl$delta_uncertainty_level))
dat_cog_ctrl$age_in_years_non_Z <- as.numeric(dat_cog_ctrl$age_in_years)
dat_cog_ctrl$age_in_years <- scale(as.numeric(dat_cog_ctrl$age_in_years))
dat_cog_ctrl$age_group_coded <- as.factor(dat_cog_ctrl$age_group_coded)
dat_cog_ctrl$RT_info_choice <- scale(as.numeric(dat_cog_ctrl$RT_info_choice))
dat_cog_ctrl$subject_ID <- as.factor(dat_cog_ctrl$subject_ID)
dat_cog_ctrl$gender_coded <- as.factor(dat_cog_ctrl$gender_coded)
contrasts(dat_cog_ctrl$gender_coded) <- contr.helmert(2)
dat_cog_ctrl$info_choice <- as.factor(dat_cog_ctrl$info_choice) # 0 = left, 1 = right
dat_cog_ctrl$percent_comprehension_non_Z <- as.numeric(dat_cog_ctrl$percent_comprehension)
dat_cog_ctrl$percent_comprehension <- scale(as.numeric(dat_cog_ctrl$percent_comprehension))
dat_cog_ctrl$chance <- 0.5

## SE2 Compute scaled deltas
# Find min and max and rescale between -1 and 1
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
dat_cog_ctrl <- rescale_deltas(dat_cog_ctrl, c("uncertainty_level"))

# Check rescaled and raw deltas have the same sign
check_delta_signs <- function(dat, deltas=c("EV", "uncertainty_level", "agency")) {
  num_trials <- length(dat$subject_ID)
  
  if ("EV" %in% deltas){
    check_EV_sign <- num_trials
    for (i in 1:num_trials) {
      if(((dat$delta_EV_non_Z[i] < 0) & (dat$delta_EV[i] < 0))
         | ((dat$delta_EV_non_Z[i] > 0) & (dat$delta_EV[i] > 0))
         | ((dat$delta_EV_non_Z[i] == 0) & (dat$delta_EV[i] == 0))) {
        check_EV_sign = check_EV_sign - 1
      }
    }
    if (check_EV_sign == 0) {
      print("rescale delta EV check passed")
    } else {
      print("rescale delta EV check failed")
    }
  }
  
  if ("uncertainty_level" %in% deltas){
    check_uncertainty_level_sign <- num_trials
    for (i in 1:num_trials) {
      if(((dat$delta_uncertainty_level_non_Z[i] < 0) & (dat$delta_uncertainty_level[i] < 0))
         | ((dat$delta_uncertainty_level_non_Z[i] > 0) & (dat$delta_uncertainty_level[i] > 0))
         | ((dat$delta_uncertainty_level_non_Z[i] == 0) & (dat$delta_uncertainty_level[i] == 0))) {
        check_uncertainty_level_sign = check_uncertainty_level_sign - 1
      }
    }
    if (check_uncertainty_level_sign == 0) {
      print("rescale delta uncertainty_level check passed")
    } else {
      print("rescale delta uncertainty_level check failed")
    }
  }
  
  if ("agency" %in% deltas) {
    check_agency_sign <- num_trials
    for (i in 1:num_trials) {
      if(((dat$delta_agency_non_Z[i] < 0) & (dat$delta_agency[i] < 0))
         | ((dat$delta_agency_non_Z[i] > 0) & (dat$delta_agency[i] > 0))
         | ((dat$delta_agency_non_Z[i] == 0) & (dat$delta_agency[i] == 0))) {
        check_agency_sign = check_agency_sign - 1
      }
    }
    if (check_agency_sign == 0) {
      print("rescale delta agency check passed")
    } else {
      print("rescale delta agency check failed")
    }
  }
}
check_delta_signs(dat_cog_ctrl, c("uncertainty_level"))

#### SE2 Demographics ####
summary(subset(dat_cog_ctrl, condition=="cognitive_1")$age_in_years_non_Z)
sdamr::sample_sd(subset(dat_cog_ctrl, condition=="cognitive_1")$age_in_years_non_Z)
table(subset(dat_cog_ctrl, condition=="cognitive_1")$gender)

#### SE2 Comprehension scores ####
mean(subset(dat_cog_ctrl, condition == "cognitive_1")$percent_comprehension_non_Z)
sdamr::sample_sd(subset(dat_cog_ctrl, condition == "cognitive_1")$percent_comprehension_non_Z)

#### SE2 Models and significance tests ####
## Models
# Full
child_cog_ctrl_mod_full <- glmer(info_choice ~ delta_uncertainty_level + 
                                   (delta_uncertainty_level
                                    | subject_ID), 
                                 data = dat_cog_ctrl, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
summary(child_cog_ctrl_mod_full)
# Without uncertainty
child_cog_ctrl_mod_full_drop_C_fixed <- glmer(info_choice ~  1 +
                                                (delta_uncertainty_level
                                                 | subject_ID), 
                                              data = dat_cog_ctrl, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
## Significance tests
child_test_fixed_uncertainty <- anova(child_cog_ctrl_mod_full, child_cog_ctrl_mod_full_drop_C_fixed, test="Chisq")
