# AUTHOR: GAIA MOLINARO
# Data analysis for the main and supplementary experiments
# in Molinaro, Cogliati Dezza, & Sharot (in prep.)

# CLEAR WORKSPACE
rm(list = ls())

# SET WORKING DIRECTORY
wd = "insert_your_wd_here"
setwd(wd)

packages <- c("ggplot2", "GGally", "ggpubr", "psych", "nls2", "RColorBrewer",
              "grid", "optimx", "ggcorrplot", "tidyverse", "glue", "gridExtra")
pacman::p_load(packages, character.only = TRUE)

# Additional basic functions
# Function to get the standard error from the mean
sem <- function(x) {
  return (sdamr::sample_sd(x)/sqrt(length(x)))
}
# Argmin function for plots
argmin <- function(x) {
  min <- min(x)
  for (i in 1:length(x)) {
    if (x[i] == min) {
      return (i)
    }
  }
}
# Function to add partial regression line to a plot
# takes a dataset and a formula, finds the intercept and coefficient
# of the first predictor, and returns a line
partial_reg_line <- function(data, f, color="black") {
  reg <- lm(data=data, formula=as.formula(f)) 
  intercept <- summary(reg)$coefficients[1, 1]
  slope <- summary(reg)$coefficients[2, 1]
  cat("formula: ", f, 
      "slope: ", round(slope, 2), 
      "p-value: ", round(summary(reg)$coefficients[2, 4], 3))
  return (geom_abline(intercept=intercept, slope=slope, color=color))
}

#### MAIN EXPERIMENT ####
#### Load and prep datasets #### 
rescale_variables <- TRUE  # whether to rescale the numeric variables
center_variables <- TRUE # whether to center numeric variables (except deltas)

# Whether to use covariates in the model
use_covariates <- FALSE

## Load children's data
dat_compt <- read.csv("data_main_experiment_children.csv")
dat_compt <- subset(dat_compt, !(condition %in% c("catch_1", "catch_2")))
dat_compt <- subset(dat_compt, !(age_in_years == 13))
dat_compt <- subset(dat_compt, catch_trials_score == 100)
dat_compt <- subset(dat_compt, wob > 0.5)

# Remove subjects for which individual fits did not converge
#dat_compt <- subset(dat_compt, ! dat_compt$subject_ID %in% c(3668306, 3898310, 3977620, 3981226))

# Remove subjects whose uncertainty coefficient was very low
#dat_compt <- subset(dat_compt, ! dat_compt$subject_ID %in% c(3802897, 3827089, 3644967))

# Create age groups
dat_compt$age_group <- "None"
dat_compt$age_group[dat_compt$age_in_years  %in% c(4, 5)] <- "4-5"
dat_compt$age_group[dat_compt$age_in_years  %in% c(6, 7)] <- "6-7"
dat_compt$age_group[dat_compt$age_in_years  %in% c(8, 9)] <- "8-9"
dat_compt$age_group[dat_compt$age_in_years  %in% c(10, 11, 12)] <- "10-12"
dat_compt$age_group <- factor(dat_compt$age_group, levels=c("4-5", "6-7", "8-9","10-12"))

if (rescale_variables) {
  # Turn variables into appropriate data types
  # Use scale() for standardized values
  dat_compt$delta_EV_non_Z <- as.numeric(dat_compt$delta_EV)
  dat_compt$delta_uncertainty_level_non_Z <- as.numeric(dat_compt$delta_uncertainty_level)
  dat_compt$delta_agency_non_Z <- as.numeric(dat_compt$delta_agency)
  dat_compt$delta_SD <- as.numeric(dat_compt$delta_SD)
  dat_compt$age_in_years_non_Z <- as.numeric(dat_compt$age_in_years)
  dat_compt$age_in_years <- scale(as.numeric(dat_compt$age_in_years), center=center_variables)
  dat_compt$age_group_coded <- as.factor(dat_compt$age_group_coded)
  dat_compt$RT_info_choice <- scale(as.numeric(dat_compt$RT_info_choice), center=center_variables)
  dat_compt$subject_ID <- as.factor(dat_compt$subject_ID)
  dat_compt$gender_coded[dat_compt$gender_coded==3] <- 1 # 1 = male, 2 = female, 3 = other
  dat_compt$gender_coded <- as.factor(dat_compt$gender_coded)
  contrasts(dat_compt$gender_coded) <- contr.helmert(2)
  dat_compt$info_choice <- as.factor(dat_compt$info_choice) # 0 = left, 1 = right
  dat_compt$percent_comprehension_non_Z <- as.numeric(dat_compt$percent_comprehension)
  dat_compt$percent_comprehension <- scale(as.numeric(dat_compt$percent_comprehension), center=center_variables)
  dat_compt$wob_non_Z <- as.numeric(dat_compt$wob)
  dat_compt$wob <- scale(as.numeric(dat_compt$wob), center=center_variables)
  dat_compt$chance <- 0.5
  # correct fishing choices based on what participants see
  dat_compt$fishing_non_Z <- as.numeric(dat_compt$prop_correct_fishing)
  dat_compt$fishing <- scale(as.numeric(dat_compt$prop_correct_fishing), center=center_variables)
} else {
  # No rescaling
  dat_compt$delta_EV_non_Z <- as.numeric(dat_compt$delta_EV)
  dat_compt$delta_uncertainty_level_non_Z <- as.numeric(dat_compt$delta_uncertainty_level)
  dat_compt$delta_agency_non_Z <- as.numeric(dat_compt$delta_agency)
  dat_compt$delta_SD <- as.numeric(dat_compt$delta_SD)
  dat_compt$age_in_years_non_Z <- as.numeric(dat_compt$age_in_years)
  dat_compt$age_in_years <- (as.numeric(dat_compt$age_in_years))
  dat_compt$age_group_coded <- as.factor(dat_compt$age_group_coded)
  dat_compt$RT_info_choice <- (as.numeric(dat_compt$RT_info_choice))
  dat_compt$subject_ID <- as.factor(dat_compt$subject_ID)
  dat_compt$gender_coded[dat_compt$gender_coded==3] <- 1 # 1 = male, 2 = female, 3 = other
  dat_compt$gender_coded <- as.factor(dat_compt$gender_coded)
  contrasts(dat_compt$gender_coded) <- contr.helmert(2)
  dat_compt$info_choice <- as.factor(dat_compt$info_choice) # 0 = left, 1 = right
  dat_compt$percent_comprehension_non_Z <- as.numeric(dat_compt$percent_comprehension)
  dat_compt$percent_comprehension <- (as.numeric(dat_compt$percent_comprehension))
  dat_compt$wob_non_Z <- as.numeric(dat_compt$wob)
  dat_compt$wob <- (as.numeric(dat_compt$wob))
  dat_compt$chance <- 0.5
  # correct fishing choices based on what participants see
  dat_compt$fishing_non_Z <- as.numeric(dat_compt$prop_correct_fishing)
  dat_compt$fishing <- (as.numeric(dat_compt$prop_correct_fishing))
}


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

# 
if (rescale_variables) {
  dat_compt <- rescale_deltas(dat_compt)  
}

# 
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

check_delta_signs(dat_compt)

# Create subsets of data for each age group
dat_compt_4_5 <- subset(dat_compt, age_group == "4-5")
dat_compt_6_7 <- subset(dat_compt, age_group == "6-7")
dat_compt_8_9 <- subset(dat_compt, age_group == "8-9")
dat_compt_10_12 <- subset(dat_compt, age_group == "10-12")
# Re-Z-score and rescale variables
if (rescale_variables) {
  dat_compt_4_5$percent_comprehension <- scale(dat_compt_4_5$percent_comprehension_non_Z, center=center_variables)
  dat_compt_6_7$percent_comprehension <- scale(dat_compt_6_7$percent_comprehension_non_Z, center=center_variables)
  dat_compt_8_9$percent_comprehension <- scale(dat_compt_8_9$percent_comprehension_non_Z, center=center_variables)
  dat_compt_10_12$percent_comprehension <- scale(dat_compt_10_12$percent_comprehension_non_Z, center=center_variables)
  dat_compt_4_5$wob <- scale(dat_compt_4_5$wob_non_Z, center=center_variables)
  dat_compt_6_7$wob <- scale(dat_compt_6_7$wob_non_Z, center=center_variables)
  dat_compt_8_9$wob <- scale(dat_compt_8_9$wob_non_Z, center=center_variables)
  dat_compt_10_12$wob <- scale(dat_compt_10_12$wob_non_Z, center=center_variables)
  dat_compt_4_5 <- rescale_deltas(dat_compt_4_5)
  dat_compt_6_7 <- rescale_deltas(dat_compt_6_7)
  dat_compt_10_12 <- rescale_deltas(dat_compt_10_12) 
}

## Load adults' data
setwd(wd)
adu_dat_compt <- read.csv("data_main_experiment_adults.csv")
adu_dat_compt <- subset(adu_dat_compt, !(condition %in% c("catch_1", "catch_2")))
adu_dat_compt <- subset(adu_dat_compt, (catch_trials_score == 100) & (wob > 0.5))

# Fix trial number
adu_dat_compt$trial_number <- substr(adu_dat_compt$condition, 7, 8)

# Turn variables into appropriate data types
# With rescaling
if (rescale_variables) {
  adu_dat_compt$delta_EV_non_Z <- as.numeric(adu_dat_compt$delta_EV)
  adu_dat_compt$delta_uncertainty_level_non_Z <- as.numeric(adu_dat_compt$delta_uncertainty_level)
  adu_dat_compt$delta_agency_non_Z <- as.numeric(adu_dat_compt$delta_agency)
  adu_dat_compt$delta_SD <- as.numeric(adu_dat_compt$delta_SD)
  adu_dat_compt$age_in_years_non_Z <- as.numeric(adu_dat_compt$age_in_years)
  adu_dat_compt$age_in_years <- scale(as.numeric(adu_dat_compt$age_in_years), center=center_variables)
  adu_dat_compt$age_group <- "Adult"
  adu_dat_compt$age_group_coded <- as.factor(4)
  adu_dat_compt$RT_info_choice <- scale(as.numeric(adu_dat_compt$RT_info_choice))
  adu_dat_compt$subject_ID <- as.factor(adu_dat_compt$subject_ID)
  adu_dat_compt$gender_coded[adu_dat_compt$gender_coded==3] <- 1 # contrast will be female vs non-female
  adu_dat_compt$gender_coded <- as.factor(adu_dat_compt$gender_coded)
  contrasts(adu_dat_compt$gender_coded) <- contr.helmert(2)
  adu_dat_compt$info_choice <- as.factor(adu_dat_compt$info_choice) # 0 = left, 1 = right
  adu_dat_compt$percent_comprehension_non_Z <- as.numeric(adu_dat_compt$percent_comprehension)
  adu_dat_compt$percent_comprehension <- scale(as.numeric(adu_dat_compt$percent_comprehension), center=center_variables)
  adu_dat_compt$wob_non_Z <- as.numeric(adu_dat_compt$wob)
  adu_dat_compt$wob <- scale(as.numeric(adu_dat_compt$wob), center=center_variables)
  adu_dat_compt$chance <- 0.5
  # correct fishing choices based on what participants see
  adu_dat_compt$fishing_non_Z <- as.numeric(adu_dat_compt$prop_correct_fishing)
  adu_dat_compt$fishing <- scale(as.numeric(adu_dat_compt$prop_correct_fishing), center=center_variables)
} else {
  adu_dat_compt$delta_EV_non_Z <- as.numeric(adu_dat_compt$delta_EV)
  adu_dat_compt$delta_uncertainty_level_non_Z <- as.numeric(adu_dat_compt$delta_uncertainty_level)
  adu_dat_compt$delta_agency_non_Z <- as.numeric(adu_dat_compt$delta_agency)
  adu_dat_compt$delta_SD <- as.numeric(adu_dat_compt$delta_SD)
  adu_dat_compt$age_in_years_non_Z <- as.numeric(adu_dat_compt$age_in_years)
  adu_dat_compt$age_in_years <- (as.numeric(adu_dat_compt$age_in_years))
  adu_dat_compt$age_group <- "Adult"
  adu_dat_compt$age_group_coded <- as.factor(4)
  adu_dat_compt$RT_info_choice <- (as.numeric(adu_dat_compt$RT_info_choice))
  adu_dat_compt$subject_ID <- as.factor(adu_dat_compt$subject_ID)
  adu_dat_compt$gender_coded[adu_dat_compt$gender_coded==3] <- 1 # contrast will be female vs non-female
  adu_dat_compt$gender_coded <- as.factor(adu_dat_compt$gender_coded)
  contrasts(adu_dat_compt$gender_coded) <- contr.helmert(2)
  adu_dat_compt$info_choice <- as.factor(adu_dat_compt$info_choice) # 0 = left, 1 = right
  adu_dat_compt$percent_comprehension_non_Z <- as.numeric(adu_dat_compt$percent_comprehension)
  adu_dat_compt$percent_comprehension <- (as.numeric(adu_dat_compt$percent_comprehension))
  adu_dat_compt$wob_non_Z <- as.numeric(adu_dat_compt$wob)
  adu_dat_compt$wob <- (as.numeric(adu_dat_compt$wob))
  adu_dat_compt$chance <- 0.5
  # correct fishing choices based on what participants see
  adu_dat_compt$fishing_non_Z <- as.numeric(adu_dat_compt$prop_correct_fishing)
  adu_dat_compt$fishing <- (as.numeric(adu_dat_compt$prop_correct_fishing))
}

## Compute scaled deltas
# Find min and max and rescale between -1 and 1
if (rescale_variables) {adu_dat_compt <- rescale_deltas(adu_dat_compt)}

# Check rescaled and raw deltas have the same sign (adults)
check_delta_signs(adu_dat_compt)

## Merged datasets (children + adults)
# Create data subsets
cols <- c("subject_ID", "info_choice", 
          "delta_EV", "delta_uncertainty_level", "delta_agency", 
          "delta_EV_non_Z", "delta_uncertainty_level_non_Z", "delta_agency_non_Z", 
          "EV_L", "EV_R", "uncertainty_level_L", "uncertainty_level_R",
          "agency_probL", "agency_probR",
          "gender_coded", "wob", "percent_comprehension", 
          "wob_non_Z", "percent_comprehension_non_Z", 
          "age_group", "age_group_coded",
          "fishing_non_Z")
dat_compt_child_adu <- rbind(dplyr::select(dat_compt, all_of(cols)), 
                             dplyr::select(adu_dat_compt, all_of(cols)))
dat_compt_4_5_adu <- rbind(dplyr::select(dat_compt_4_5, all_of(cols)), 
                           dplyr::select(adu_dat_compt, all_of(cols)))
dat_compt_6_7_adu <- rbind(dplyr::select(dat_compt_6_7, all_of(cols)), 
                           dplyr::select(adu_dat_compt, all_of(cols)))
dat_compt_8_9_adu <- rbind(dplyr::select(dat_compt_8_9, all_of(cols)), 
                           dplyr::select(adu_dat_compt, all_of(cols)))
dat_compt_10_12_adu <- rbind(dplyr::select(dat_compt_10_12, all_of(cols)), 
                             dplyr::select(adu_dat_compt, all_of(cols)))


# Rescale variables
if (rescale_variables) {
  dat_compt_child_adu$percent_comprehension <- scale(dat_compt_child_adu$percent_comprehension_non_Z, center=center_variables)
  dat_compt_4_5_adu$percent_comprehension <- scale(dat_compt_4_5_adu$percent_comprehension_non_Z, center=center_variables)
  dat_compt_6_7_adu$percent_comprehension <- scale(dat_compt_6_7_adu$percent_comprehension_non_Z, center=center_variables)
  dat_compt_8_9_adu$percent_comprehension <- scale(dat_compt_8_9_adu$percent_comprehension_non_Z, center=center_variables)
  dat_compt_10_12_adu$percent_comprehension <- scale(dat_compt_10_12_adu$percent_comprehension_non_Z, center=center_variables)
  dat_compt_child_adu$wob <- scale(dat_compt_child_adu$wob_non_Z, center=center_variables)
  dat_compt_4_5_adu$wob <- scale(dat_compt_4_5_adu$wob_non_Z, center=center_variables)
  dat_compt_6_7_adu$wob <- scale(dat_compt_6_7_adu$wob_non_Z, center=center_variables)
  dat_compt_8_9_adu$wob <- scale(dat_compt_8_9_adu$wob_non_Z, center=center_variables)
  dat_compt_10_12_adu$wob <- scale(dat_compt_10_12_adu$wob_non_Z, center=center_variables)
  dat_compt_child_adu <- rescale_deltas(dat_compt_child_adu)
  dat_compt_4_5_adu <- rescale_deltas(dat_compt_4_5_adu)
  dat_compt_6_7_adu <- rescale_deltas(dat_compt_6_7_adu)
  dat_compt_8_9_adu <- rescale_deltas(dat_compt_8_9_adu)
  dat_compt_10_12_adu <- rescale_deltas(dat_compt_10_12_adu) 
}

# Set gender contrasts
contrasts(dat_compt_child_adu$gender_coded) <- contr.helmert(2)
contrasts(dat_compt_4_5_adu$gender_coded) <- contr.helmert(2)
contrasts(dat_compt_6_7_adu$gender_coded) <- contr.helmert(2)
contrasts(dat_compt_8_9_adu$gender_coded) <- contr.helmert(2)
contrasts(dat_compt_10_12_adu$gender_coded) <- contr.helmert(2)

# Create new variable and  contrasts for children vs adults factor
dat_compt_child_adu$group <- ifelse(dat_compt_child_adu$age_group_coded == 4, 1, -1) # -1 = child, 1 = adult
dat_compt_4_5_adu$group <- ifelse(dat_compt_4_5_adu$age_group_coded == 4,  1, -1) # -1 = child, 1 = adult
dat_compt_6_7_adu$group <-  ifelse(dat_compt_6_7_adu$age_group_coded == 4, 1, -1) # -1 = child, 1 = adult
dat_compt_8_9_adu$group <-  ifelse(dat_compt_8_9_adu$age_group_coded == 4, 1, -1) # -1 = child, 1 = adult
dat_compt_10_12_adu$group <-  ifelse(dat_compt_10_12_adu$age_group_coded == 4, 1, -1) # -1 = child, 1 = adult

#### Demographics ####
summary(subset(dat_compt, trial_number==1)$age_in_years_non_Z)
sdamr::sample_sd(subset(dat_compt, trial_number==1)$age_in_years_non_Z)
psych::describeBy(subset(dat_compt, trial_number==1)$age_in_years_non_Z, 
                  subset(dat_compt, trial_number==1)$age_group)
table(subset(dat_compt, trial_number==1)$gender)
table(subset(dat_compt, trial_number==1)$gender, 
      subset(dat_compt, trial_number==1)$age_group)
psych::describeBy(subset(adu_dat_compt, trial_number==1)$age_in_years_non_Z, 
                  subset(adu_dat_compt, trial_number==1)$age_group)
table(subset(adu_dat_compt, trial_number==1)$gender, 
      subset(adu_dat_compt, trial_number==1)$age_group)

#### Comprehension and wob scores ####
## t-tests
t.test(subset(dat_compt, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(adu_dat_compt, trial_number == 1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_compt_4_5, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_compt_6_7, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_compt_8_9, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_compt_10_12, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)

# Standard deviations
sdamr::sample_sd(subset(dat_compt, trial_number ==1)$percent_comprehension_non_Z)
sdamr::sample_sd(subset(adu_dat_compt, trial_number ==1)$percent_comprehension_non_Z)
sdamr::sample_sd(subset(dat_compt, trial_number ==1)$wob_non_Z)
sdamr::sample_sd(subset(adu_dat_compt, trial_number ==1)$wob_non_Z)

# Descriptive stats by age group
psych::describeBy(subset(dat_compt, trial_number==1)$percent_comprehension_non_Z, 
                  subset(dat_compt, trial_number==1)$age_group)
psych::describeBy(subset(dat_compt, trial_number==1)$wob_non_Z, 
                  subset(dat_compt, trial_number==1)$age_group)

## Create long form data sets
# Create separate data sets with variables for children and adults' pc and wob
dat_compt_pc_wob <- dplyr::select(subset(dat_compt, trial_number==1), 
                                  all_of(c("subject_ID", "age_group_coded", 
                                           "percent_comprehension_non_Z", 
                                           "wob_non_Z")))
dat_compt_pc_wob$percent_comprehension_child <- dat_compt_pc_wob$"percent_comprehension_non_Z"
dat_compt_pc_wob$wob_child <- dat_compt_pc_wob$"wob_non_Z"

adu_dat_compt_pc_wob <- dplyr::select(subset(adu_dat_compt, trial_number==1), 
                                      all_of(c("subject_ID", "age_group_coded", 
                                               "percent_comprehension_non_Z", 
                                               "wob_non_Z")))

# Merge datasets
dat_compt_child_adu_pc_wob <- rbind(dplyr::select(dat_compt_pc_wob, all_of(
  c("subject_ID", "age_group_coded", "percent_comprehension_non_Z", "wob_non_Z"))), 
  dplyr::select(adu_dat_compt_pc_wob, all_of(
    c("subject_ID", "age_group_coded", "percent_comprehension_non_Z", "wob_non_Z"))))

# Create group factor 
dat_compt_child_adu_pc_wob$group <- as.factor(ifelse(dat_compt_child_adu_pc_wob$age_group_coded == 4, "Adults", "Children"))

# Turn into long form 
dat_compt_child_adu_pc_wob <- tidyr::gather(dat_compt_child_adu_pc_wob, 
                                            key = "variable", value = "score", percent_comprehension_non_Z:wob_non_Z, factor_key=TRUE)
# Create a new variable to distinguish between children's and adult's scores
dat_compt_child_adu_pc_wob$variable2 <- "None"
dat_compt_child_adu_pc_wob$variable2[dat_compt_child_adu_pc_wob$group == "Children" 
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_child"
dat_compt_child_adu_pc_wob$variable2[dat_compt_child_adu_pc_wob$group == "Children" 
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_child"
dat_compt_child_adu_pc_wob$variable2[dat_compt_child_adu_pc_wob$group == "Adults" 
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_adu"
dat_compt_child_adu_pc_wob$variable2[dat_compt_child_adu_pc_wob$group == "Adults" 
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_adu"
dat_compt_child_adu_pc_wob$variable2 <- factor(dat_compt_child_adu_pc_wob$variable2, 
                                               levels=c("percent_comprehension_child",
                                                        "percent_comprehension_adu",
                                                        "wob_child",
                                                        "wob_adu"))
# Create a new variable for colors for variable2
dat_compt_child_adu_pc_wob$color <- "None"
dat_compt_child_adu_pc_wob$color <- as.factor(ifelse(
  dat_compt_child_adu_pc_wob$age_group_coded == 4, "#444444", "#AAAAAA"))

# Create a new variable to distinguish between each group of children's and adult's scores
dat_compt_child_adu_pc_wob$variable3 <- "None"
dat_compt_child_adu_pc_wob$variable3[dat_compt_child_adu_pc_wob$age_group_coded == 0
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_4_5"
dat_compt_child_adu_pc_wob$variable3[dat_compt_child_adu_pc_wob$age_group_coded == 0
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_4_5"
dat_compt_child_adu_pc_wob$variable3[dat_compt_child_adu_pc_wob$age_group_coded == 1
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_6_7"
dat_compt_child_adu_pc_wob$variable3[dat_compt_child_adu_pc_wob$age_group_coded == 1
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_6_7"
dat_compt_child_adu_pc_wob$variable3[dat_compt_child_adu_pc_wob$age_group_coded == 2
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_8_9"
dat_compt_child_adu_pc_wob$variable3[dat_compt_child_adu_pc_wob$age_group_coded == 2
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_8_9"
dat_compt_child_adu_pc_wob$variable3[dat_compt_child_adu_pc_wob$age_group_coded == 3
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_10_12"
dat_compt_child_adu_pc_wob$variable3[dat_compt_child_adu_pc_wob$age_group_coded == 3
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_10_12"
dat_compt_child_adu_pc_wob$variable3[dat_compt_child_adu_pc_wob$age_group_coded == 4 
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_adu"
dat_compt_child_adu_pc_wob$variable3[dat_compt_child_adu_pc_wob$age_group_coded == 4
                                     & dat_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_adu"
dat_compt_child_adu_pc_wob$variable3 <- factor(dat_compt_child_adu_pc_wob$variable3, 
                                               levels=c("percent_comprehension_4_5",
                                                        "percent_comprehension_6_7",
                                                        "percent_comprehension_8_9",
                                                        "percent_comprehension_10_12",
                                                        "percent_comprehension_adu",
                                                        "wob_4_5",
                                                        "wob_6_7",
                                                        "wob_8_9",
                                                        "wob_10_12",
                                                        "wob_adu"))

table(dat_compt_child_adu_pc_wob$variable3)

## Plots
# By group
compt_pc_by_group <- ggplot(subset(dat_compt_child_adu_pc_wob, variable2 %in% 
                                     c("percent_comprehension_child", "percent_comprehension_adu")), 
                            aes(x=variable3, y=score, fill=variable2)) +
  geom_jitter(height=0, color = subset(dat_compt_child_adu_pc_wob, 
                                       variable2 %in% c("percent_comprehension_child", "percent_comprehension_adu"))$color,
              size=1.5, alpha=0.9)+   
  geom_violin(alpha=0.3) +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0, 1)) + 
  scale_fill_manual(values=c("#AAAAAA", "#444444")) +
  #stat_summary(fun=mean, geom="point", shape=20, size=4, color="black", fill="black") +
  theme_classic() + labs(x="comprehension", y="Score") + 
  theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm")) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

compt_wob_by_group <- ggplot(subset(dat_compt_child_adu_pc_wob, variable2 %in% 
                                      c("wob_child", "wob_adu")), 
                             aes(x=variable3, y=score, fill=variable2)) +
  geom_jitter(height=0, color = subset(dat_compt_child_adu_pc_wob, 
                                       variable2 %in% c("wob_child", "wob_adu"))$color,
              size=1.5, alpha=0.9)+   
  geom_violin(alpha=0.3) +
  geom_hline(yintercept=0.5, color="#AAAAAA") +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0, 1)) + 
  scale_fill_manual(values=c("#AAAAAA", "#444444")) +
  #stat_summary(fun=mean, geom="point", shape=20, size=4, color="black", fill="black") +
  theme_classic() + labs(x="comprehension", y="Score") + 
  theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm")) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

grid.arrange(compt_pc_by_group, compt_wob_by_group, ncol=2)


#### Fishing choices ####
# t-tests

# Overall fishing choices compared to chance and children vs adults
# based on what participants see
sdamr::sample_sd(subset(dat_compt, trial_number == 1)$fishing_non_Z)
t.test(subset(adu_dat_compt, trial_number == 1)$fishing_non_Z, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(adu_dat_compt, trial_number == 1)$fishing_non_Z)
t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_non_Z, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_non_Z)
t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_non_Z, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_non_Z)
t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_non_Z, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_non_Z)
t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_non_Z, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_non_Z)

# Plot average scores by age group
# Create dataset for fishing choices
age_group <- factor(c("4-5", "6-7", "8-9", "10-12", "Adults"), levels=c("4-5", "6-7", "8-9", "10-12", "Adults"))
x_points <- c(0.1, 0.2, 0.3, 0.4, 0.5)
x_labels <- c("4-5", "6-7", "8-9", "10-12", "Adults")

# based on what participants see
fishing <- c(
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_non_Z),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_non_Z),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_non_Z),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_non_Z),
  mean(subset(adu_dat_compt, trial_number == 1)$fishing_non_Z)
)
fishing_sem <- c(
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_non_Z),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_non_Z),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_non_Z),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_non_Z),
  sem(subset(adu_dat_compt, trial_number == 1)$fishing_non_Z)
)

fishing_choices_dat_by_group <- data.frame(x_points, x_labels,
                                           fishing, fishing_sem)

fishing_choices_by_group_plot <- ggplot(fishing_choices_dat_by_group) + 
  geom_line(aes(y=fishing, x=x_points), size=1, color="darkgrey")+
  geom_point(aes(y=fishing, x=x_points), color="darkgrey", size = 4)+
  geom_ribbon(aes(ymin=fishing-fishing_sem, 
                  ymax=fishing+fishing_sem, 
                  x=x_points), fill="darkgrey" , alpha = 0.3)+
  theme_classic()+
  theme(text = element_text(size=16), legend.position="none") +
  scale_x_continuous(labels=fishing_choices_dat_by_group$x_labels) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  geom_hline(yintercept=0.5) +
  labs(x="Age group", y="Correct fishing choices\n(based on what participants see)")


#### Correlations between deltas ####
corr_variables <- c("delta_EV", "delta_SD", "delta_uncertainty_level", "delta_agency")
dat_compt_corr <- dat_compt[corr_variables]
dat_compt_corr_mat <- round(cor(dat_compt_corr, method="pearson"), 3)
dat_compt_p_mat <- ggcorrplot::cor_pmat(dat_compt_corr)
ggcorrplot(dat_compt_corr_mat, outline.col = "white", lab = TRUE, insig = "blank")
dat_compt_corr_mat 
dat_compt_p_mat

adu_dat_compt_corr <- adu_dat_compt[corr_variables]
adu_dat_compt_corr_mat <- round(cor(adu_dat_compt_corr, method="pearson"), 3)
adu_dat_compt_p_mat <- ggcorrplot::cor_pmat(adu_dat_compt_corr)
ggcorrplot(adu_dat_compt_corr_mat, outline.col = "white", lab = TRUE, insig = "blank")
adu_dat_compt_corr_mat 
adu_dat_compt_p_mat

####  Models by age group #### 
if (use_covariates) {
  # USING COVARIATES
  # Adults
  # Chance only
  adu_compt_mod_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                                data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full 
  adu_compt_mod_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + gender_coded + 
                                (delta_EV + delta_uncertainty_level + delta_agency 
                                 | subject_ID), 
                              data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without intercept
  adu_compt_mod_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                     percent_comprehension + wob + gender_coded -1 +
                                                     (delta_EV + delta_uncertainty_level + delta_agency 
                                                      | subject_ID), 
                                                   data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without EV as a fixed effect
  adu_compt_mod_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  adu_compt_mod_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                             percent_comprehension + wob + gender_coded + 
                                             (delta_EV + delta_uncertainty_level + delta_agency 
                                              | subject_ID), 
                                           data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  adu_compt_mod_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                              percent_comprehension + wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without comprehension as a fixed effect
  adu_compt_mod_full_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              wob + gender_coded +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without wob as a fixed effect
  adu_compt_mod_full_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + gender_coded +
                                               (delta_EV + delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without gender as a fixed effect
  adu_compt_mod_full_drop_gender_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                  percent_comprehension + wob  +
                                                  (delta_EV + delta_uncertainty_level + delta_agency
                                                   | subject_ID), 
                                                data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop EV random
  adu_compt_mod_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob  + gender_coded +
                                               (delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  adu_compt_mod_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob  + gender_coded +
                                              (delta_EV + delta_agency
                                               | subject_ID), 
                                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  adu_compt_mod_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob  + gender_coded +
                                               (delta_EV + delta_uncertainty_level
                                                | subject_ID), 
                                             data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Affect only
  adu_compt_mod_Af <- glmer(info_choice ~ delta_EV +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV 
                               | subject_ID), 
                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only
  adu_compt_mod_C <- glmer(info_choice ~  delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + 
                             (delta_uncertainty_level 
                              | subject_ID), 
                           data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only 
  adu_compt_mod_Ac <- glmer(info_choice ~ delta_agency +
                              percent_comprehension + wob + gender_coded +
                              (delta_agency 
                               | subject_ID), 
                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition 
  adu_compt_mod_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               percent_comprehension + wob + gender_coded + 
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID), 
                             data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  adu_compt_mod_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                percent_comprehension + wob + gender_coded + 
                                (delta_EV + delta_agency
                                 | subject_ID), 
                              data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  adu_compt_mod_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + 
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID), 
                             data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Children's age groups
  # 4-5
  # Chance only
  compt_mod_4_5_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                                data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full
  compt_mod_4_5_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + gender_coded + 
                                (delta_EV + delta_uncertainty_level + delta_agency 
                                 | subject_ID), 
                              data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without intercept
  compt_mod_4_5_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                     percent_comprehension + wob + gender_coded -1 + 
                                                     (delta_EV + delta_uncertainty_level + delta_agency 
                                                      | subject_ID), 
                                                   data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without EV as a fixed effect
  compt_mod_4_5_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  compt_mod_4_5_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                             percent_comprehension + wob + gender_coded + 
                                             (delta_EV + delta_uncertainty_level + delta_agency 
                                              | subject_ID), 
                                           data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  compt_mod_4_5_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                              percent_comprehension + wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without pc as a fixed effect
  compt_mod_4_5_full_drop_pc_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                              wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without wob as a fixed effect
  compt_mod_4_5_full_drop_wob_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + gender_coded + 
                                               (delta_EV + delta_uncertainty_level + delta_agency 
                                                | subject_ID), 
                                             data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without gender as a fixed effect
  compt_mod_4_5_full_drop_gender_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                                  percent_comprehension +  wob + 
                                                  (delta_EV + delta_uncertainty_level + delta_agency 
                                                   | subject_ID), 
                                                data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect only
  compt_mod_4_5_Af <- glmer(info_choice ~ delta_EV +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV 
                               | subject_ID), 
                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only
  compt_mod_4_5_C <- glmer(info_choice ~  delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + 
                             (delta_uncertainty_level 
                              | subject_ID), 
                           data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only
  compt_mod_4_5_Ac <- glmer(info_choice ~ delta_agency +
                              percent_comprehension + wob + gender_coded +
                              (delta_agency 
                               | subject_ID), 
                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition
  compt_mod_4_5_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               percent_comprehension + wob + gender_coded + 
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID), 
                             data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  compt_mod_4_5_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                percent_comprehension + wob + gender_coded + 
                                (delta_EV + delta_agency
                                 | subject_ID), 
                              data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  compt_mod_4_5_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + 
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID), 
                             data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Drop EV random
  compt_mod_4_5_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob  + gender_coded +
                                               (delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  compt_mod_4_5_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob  + gender_coded +
                                              (delta_EV + delta_agency
                                               | subject_ID), 
                                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  compt_mod_4_5_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob  + gender_coded +
                                               (delta_EV + delta_uncertainty_level
                                                | subject_ID), 
                                             data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  
  # 6-7
  # Chance only
  compt_mod_6_7_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                                data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full 
  compt_mod_6_7_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + gender_coded + 
                                (delta_EV + delta_uncertainty_level + delta_agency
                                 | subject_ID), 
                              data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without intercept
  compt_mod_6_7_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                     percent_comprehension + wob + gender_coded -1 + 
                                                     (delta_EV + delta_uncertainty_level + delta_agency 
                                                      | subject_ID), 
                                                   data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without EV as a fixed effect
  compt_mod_6_7_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  compt_mod_6_7_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                             percent_comprehension + wob + gender_coded + 
                                             (delta_EV + delta_uncertainty_level + delta_agency 
                                              | subject_ID), 
                                           data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  compt_mod_6_7_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                              percent_comprehension + wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without pc as a fixed effect
  compt_mod_6_7_full_drop_pc_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                              wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without wob as a fixed effect
  compt_mod_6_7_full_drop_wob_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + gender_coded + 
                                               (delta_EV + delta_uncertainty_level + delta_agency 
                                                | subject_ID), 
                                             data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without gender as a fixed effect
  compt_mod_6_7_full_drop_gender_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                                  percent_comprehension +  wob + 
                                                  (delta_EV + delta_uncertainty_level + delta_agency 
                                                   | subject_ID), 
                                                data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Affect only
  compt_mod_6_7_Af <- glmer(info_choice ~ delta_EV +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV
                               | subject_ID), 
                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only 
  compt_mod_6_7_C <- glmer(info_choice ~  delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + 
                             (delta_uncertainty_level 
                              | subject_ID), 
                           data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only
  compt_mod_6_7_Ac <- glmer(info_choice ~ delta_agency +
                              percent_comprehension + wob + gender_coded +
                              (delta_agency 
                               | subject_ID), 
                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition 
  compt_mod_6_7_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               percent_comprehension + wob + gender_coded + 
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID), 
                             data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  compt_mod_6_7_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                percent_comprehension + wob + gender_coded + 
                                (delta_EV + delta_agency
                                 | subject_ID), 
                              data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  compt_mod_6_7_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + 
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID), 
                             data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop EV random
  compt_mod_6_7_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob  + gender_coded +
                                               (delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  compt_mod_6_7_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob  + gender_coded +
                                              (delta_EV + delta_agency
                                               | subject_ID), 
                                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  compt_mod_6_7_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob  + gender_coded +
                                               (delta_EV + delta_uncertainty_level
                                                | subject_ID), 
                                             data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # 8-9
  # Chance only
  compt_mod_8_9_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                                data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without intercept
  compt_mod_8_9_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                     percent_comprehension + wob + gender_coded -1 + 
                                                     (delta_EV + delta_uncertainty_level + delta_agency 
                                                      | subject_ID), 
                                                   data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full 
  compt_mod_8_9_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + gender_coded + 
                                (delta_EV + delta_uncertainty_level + delta_agency
                                 | subject_ID), 
                              data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without EV as a fixed effect
  compt_mod_8_9_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  compt_mod_8_9_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                             percent_comprehension + wob + gender_coded + 
                                             (delta_EV + delta_uncertainty_level + delta_agency 
                                              | subject_ID), 
                                           data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  compt_mod_8_9_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                              percent_comprehension + wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without pc as a fixed effect
  compt_mod_8_9_full_drop_pc_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                              wob + gender_coded + 
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without wob as a fixed effect
  compt_mod_8_9_full_drop_wob_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + gender_coded + 
                                               (delta_EV + delta_uncertainty_level + delta_agency 
                                                | subject_ID), 
                                             data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without gender as a fixed effect
  compt_mod_8_9_full_drop_gender_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                                  percent_comprehension +  wob + 
                                                  (delta_EV + delta_uncertainty_level + delta_agency 
                                                   | subject_ID), 
                                                data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Affect only
  compt_mod_8_9_Af <- glmer(info_choice ~ delta_EV +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV 
                               | subject_ID), 
                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only
  compt_mod_8_9_C <- glmer(info_choice ~  delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + 
                             (delta_uncertainty_level 
                              | subject_ID), 
                           data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only
  compt_mod_8_9_Ac <- glmer(info_choice ~ delta_agency +
                              percent_comprehension + wob + gender_coded +
                              (delta_agency 
                               | subject_ID), 
                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition
  compt_mod_8_9_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               percent_comprehension + wob + gender_coded + 
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID), 
                             data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  compt_mod_8_9_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                percent_comprehension + wob + gender_coded + 
                                (delta_EV + delta_agency
                                 | subject_ID), 
                              data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  compt_mod_8_9_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + 
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID), 
                             data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop EV random
  compt_mod_8_9_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob  + gender_coded +
                                               (delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  compt_mod_8_9_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob  + gender_coded +
                                              (delta_EV + delta_agency
                                               | subject_ID), 
                                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  compt_mod_8_9_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob  + gender_coded +
                                               (delta_EV + delta_uncertainty_level
                                                | subject_ID), 
                                             data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  
  # 10-12
  # Chance only
  compt_mod_10_12_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                                  data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full 
  compt_mod_10_12_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                  percent_comprehension + wob + gender_coded + 
                                  (delta_EV + delta_uncertainty_level + delta_agency
                                   | subject_ID), 
                                data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without intercept
  compt_mod_10_12_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                       percent_comprehension + wob + gender_coded -1 + 
                                                       (delta_EV + delta_uncertainty_level + delta_agency 
                                                        | subject_ID), 
                                                     data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Full, without EV as a fixed effect
  compt_mod_10_12_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                                percent_comprehension + wob + gender_coded + 
                                                (delta_EV + delta_uncertainty_level + delta_agency 
                                                 | subject_ID), 
                                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  compt_mod_10_12_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                               percent_comprehension + wob + gender_coded + 
                                               (delta_EV + delta_uncertainty_level + delta_agency 
                                                | subject_ID), 
                                             data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  compt_mod_10_12_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                                percent_comprehension + wob + gender_coded + 
                                                (delta_EV + delta_uncertainty_level + delta_agency 
                                                 | subject_ID), 
                                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without pc as a fixed effect
  compt_mod_10_12_full_drop_pc_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                                wob + gender_coded + 
                                                (delta_EV + delta_uncertainty_level + delta_agency 
                                                 | subject_ID), 
                                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without wob as a fixed effect
  compt_mod_10_12_full_drop_wob_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                                 percent_comprehension + gender_coded + 
                                                 (delta_EV + delta_uncertainty_level + delta_agency 
                                                  | subject_ID), 
                                               data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without gender as a fixed effect
  compt_mod_10_12_full_drop_gender_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level + delta_agency +
                                                    percent_comprehension +  wob + 
                                                    (delta_EV + delta_uncertainty_level + delta_agency 
                                                     | subject_ID), 
                                                  data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Affect only 
  compt_mod_10_12_Af <- glmer(info_choice ~ delta_EV +
                                percent_comprehension + wob + gender_coded + 
                                (delta_EV 
                                 | subject_ID), 
                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only
  compt_mod_10_12_C <- glmer(info_choice ~  delta_uncertainty_level +
                               percent_comprehension + wob + gender_coded + 
                               (delta_uncertainty_level 
                                | subject_ID), 
                             data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only 
  compt_mod_10_12_Ac <- glmer(info_choice ~ delta_agency +
                                percent_comprehension + wob + gender_coded +
                                (delta_agency 
                                 | subject_ID), 
                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition 
  compt_mod_10_12_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                 percent_comprehension + wob + gender_coded + 
                                 (delta_EV + delta_uncertainty_level 
                                  | subject_ID), 
                               data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  compt_mod_10_12_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                  percent_comprehension + wob + gender_coded + 
                                  (delta_EV + delta_agency
                                   | subject_ID), 
                                data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  compt_mod_10_12_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                                 percent_comprehension + wob + gender_coded + 
                                 (delta_uncertainty_level + delta_agency 
                                  | subject_ID), 
                               data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  
  # Drop EV random
  compt_mod_10_12_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                 percent_comprehension + wob  + gender_coded +
                                                 (delta_uncertainty_level + delta_agency
                                                  | subject_ID), 
                                               data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  compt_mod_10_12_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                percent_comprehension + wob  + gender_coded +
                                                (delta_EV + delta_agency
                                                 | subject_ID), 
                                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  compt_mod_10_12_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                 percent_comprehension + wob  + gender_coded +
                                                 (delta_EV + delta_uncertainty_level
                                                  | subject_ID), 
                                               data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
} else {
  # NO COVARIATES
  # Adults
  # Chance only
  adu_compt_mod_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                                data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full 
  adu_compt_mod_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                (delta_EV + delta_uncertainty_level + delta_agency 
                                 | subject_ID), 
                              data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without intercept
  adu_compt_mod_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                     (delta_EV + delta_uncertainty_level + delta_agency 
                                                      | subject_ID), 
                                                   data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without EV as a fixed effect
  adu_compt_mod_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  adu_compt_mod_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                             (delta_EV + delta_uncertainty_level + delta_agency 
                                              | subject_ID), 
                                           data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  adu_compt_mod_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop EV random
  adu_compt_mod_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               (delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  adu_compt_mod_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              (delta_EV + delta_agency
                                               | subject_ID), 
                                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  adu_compt_mod_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               (delta_EV + delta_uncertainty_level
                                                | subject_ID), 
                                             data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Affect only
  adu_compt_mod_Af <- glmer(info_choice ~ delta_EV +
                              (delta_EV 
                               | subject_ID), 
                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only
  adu_compt_mod_C <- glmer(info_choice ~  delta_uncertainty_level +
                             (delta_uncertainty_level 
                              | subject_ID), 
                           data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only 
  adu_compt_mod_Ac <- glmer(info_choice ~ delta_agency +
                              (delta_agency 
                               | subject_ID), 
                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition 
  adu_compt_mod_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID), 
                             data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  adu_compt_mod_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                (delta_EV + delta_agency
                                 | subject_ID), 
                              data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  adu_compt_mod_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID), 
                             data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Children's age groups
  # 4-5
  # Chance only
  compt_mod_4_5_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                                data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full
  compt_mod_4_5_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                              (delta_EV + delta_uncertainty_level + delta_agency
                                 | subject_ID), 
                              data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
   
  # Full, without intercept
  compt_mod_4_5_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                     (delta_EV + delta_uncertainty_level + delta_agency 
                                                      | subject_ID), 
                                                   data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without EV as a fixed effect
  compt_mod_4_5_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  compt_mod_4_5_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                             (delta_EV + delta_uncertainty_level + delta_agency 
                                              | subject_ID), 
                                           data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  compt_mod_4_5_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect only
  compt_mod_4_5_Af <- glmer(info_choice ~ delta_EV +
                              (delta_EV 
                               | subject_ID), 
                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only
  compt_mod_4_5_C <- glmer(info_choice ~  delta_uncertainty_level +
                             (delta_uncertainty_level 
                              | subject_ID), 
                           data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only
  compt_mod_4_5_Ac <- glmer(info_choice ~ delta_agency +
                              (delta_agency 
                               | subject_ID), 
                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition
  compt_mod_4_5_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID), 
                             data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  compt_mod_4_5_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                (delta_EV + delta_agency
                                 | subject_ID), 
                              data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  compt_mod_4_5_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID), 
                             data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Drop EV random
  compt_mod_4_5_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               (delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  compt_mod_4_5_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              (delta_EV + delta_agency
                                               | subject_ID), 
                                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  compt_mod_4_5_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               (delta_EV + delta_uncertainty_level
                                                | subject_ID), 
                                             data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  
  # 6-7
  # Chance only
  compt_mod_6_7_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                                data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full 
  compt_mod_6_7_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                (delta_EV + delta_uncertainty_level + delta_agency                               
                                 | subject_ID), 
                              data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without intercept
  compt_mod_6_7_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                     (delta_EV + delta_uncertainty_level + delta_agency 
                                                      | subject_ID), 
                                                   data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without EV as a fixed effect
  compt_mod_6_7_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  compt_mod_6_7_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                             (delta_EV + delta_uncertainty_level + delta_agency 
                                              | subject_ID), 
                                           data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  compt_mod_6_7_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect only
  compt_mod_6_7_Af <- glmer(info_choice ~ delta_EV +
                              (delta_EV
                               | subject_ID), 
                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only 
  compt_mod_6_7_C <- glmer(info_choice ~  delta_uncertainty_level +
                             (delta_uncertainty_level 
                              | subject_ID), 
                           data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only
  compt_mod_6_7_Ac <- glmer(info_choice ~ delta_agency +
                              (delta_agency 
                               | subject_ID), 
                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition 
  compt_mod_6_7_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID), 
                             data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  compt_mod_6_7_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                (delta_EV + delta_agency
                                 | subject_ID), 
                              data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  compt_mod_6_7_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID), 
                             data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop EV random
  compt_mod_6_7_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               (delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  compt_mod_6_7_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              (delta_EV + delta_agency
                                               | subject_ID), 
                                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  compt_mod_6_7_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               (delta_EV + delta_uncertainty_level
                                                | subject_ID), 
                                             data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # 8-9
  # Chance only
  compt_mod_8_9_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                                data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full, without intercept
  compt_mod_8_9_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                     (delta_EV + delta_uncertainty_level + delta_agency 
                                                      | subject_ID), 
                                                   data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full 
  compt_mod_8_9_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                (delta_EV + delta_uncertainty_level + delta_agency
                                 | subject_ID), 
                              data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without EV as a fixed effect
  compt_mod_8_9_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  compt_mod_8_9_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                             (delta_EV + delta_uncertainty_level + delta_agency 
                                              | subject_ID), 
                                           data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  compt_mod_8_9_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                              (delta_EV + delta_uncertainty_level + delta_agency 
                                               | subject_ID), 
                                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Affect only
  compt_mod_8_9_Af <- glmer(info_choice ~ delta_EV +
                              (delta_EV 
                               | subject_ID), 
                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only
  compt_mod_8_9_C <- glmer(info_choice ~  delta_uncertainty_level +
                             (delta_uncertainty_level 
                              | subject_ID), 
                           data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only
  compt_mod_8_9_Ac <- glmer(info_choice ~ delta_agency +
                              (delta_agency 
                               | subject_ID), 
                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition
  compt_mod_8_9_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID), 
                             data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  compt_mod_8_9_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                (delta_EV + delta_agency
                                 | subject_ID), 
                              data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  compt_mod_8_9_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID), 
                             data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop EV random
  compt_mod_8_9_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               (delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  compt_mod_8_9_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              (delta_EV + delta_agency
                                               | subject_ID), 
                                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  compt_mod_8_9_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               (delta_EV + delta_uncertainty_level
                                                | subject_ID), 
                                             data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  
  # 10-12
  # Chance only
  compt_mod_10_12_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                                  data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Full 
  compt_mod_10_12_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                  (delta_EV + delta_uncertainty_level + delta_agency
                                   | subject_ID), 
                                data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without intercept
  compt_mod_10_12_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                       (delta_EV + delta_uncertainty_level + delta_agency 
                                                        | subject_ID), 
                                                     data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  # Full, without EV as a fixed effect
  compt_mod_10_12_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                                (delta_EV + delta_uncertainty_level + delta_agency 
                                                 | subject_ID), 
                                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without uncertainty as a fixed effect
  compt_mod_10_12_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                               (delta_EV + delta_uncertainty_level + delta_agency 
                                                | subject_ID), 
                                             data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  # Full, without agency as a fixed effect
  compt_mod_10_12_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                                (delta_EV + delta_uncertainty_level + delta_agency 
                                                 | subject_ID), 
                                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect only 
  compt_mod_10_12_Af <- glmer(info_choice ~ delta_EV +
                                (delta_EV 
                                 | subject_ID), 
                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition only
  compt_mod_10_12_C <- glmer(info_choice ~  delta_uncertainty_level +
                               (delta_uncertainty_level 
                                | subject_ID), 
                             data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Action only 
  compt_mod_10_12_Ac <- glmer(info_choice ~ delta_agency +
                                (delta_agency 
                                 | subject_ID), 
                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Cognition 
  compt_mod_10_12_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                 (delta_EV + delta_uncertainty_level 
                                  | subject_ID), 
                               data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Affect + Action 
  compt_mod_10_12_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                  (delta_EV + delta_agency
                                   | subject_ID), 
                                data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Cognition + Action 
  compt_mod_10_12_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                                 (delta_uncertainty_level + delta_agency 
                                  | subject_ID), 
                               data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  
  # Drop EV random
  compt_mod_10_12_full_drop_Af_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                 (delta_uncertainty_level + delta_agency
                                                  | subject_ID), 
                                               data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop uncertainty level random
  compt_mod_10_12_full_drop_C_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                (delta_EV + delta_agency
                                                 | subject_ID), 
                                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  # Drop agency random
  compt_mod_10_12_full_drop_Ac_random <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                 (delta_EV + delta_uncertainty_level
                                                  | subject_ID), 
                                               data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
} # End of if/else statement for using covariates or not  

#### Significance tests within groups ####
# (within single groups)
# access results of test a with a$Chisq, a$Df, a$Pr[2]

## Adults overall
adu_test_fixed_EV <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_Af_fixed, test="Chisq")
adu_test_fixed_uncertainty <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_C_fixed, test="Chisq")
adu_test_fixed_agency <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_Ac_fixed, test="Chisq")
# covariates
adu_test_fixed_pc <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_pc_fixed, test="Chisq")
adu_test_fixed_wob <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_wob_fixed, test="Chisq")
adu_test_fixed_gender <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_gender_fixed, test="Chisq")

## 4-5
test_fixed_4_5_EV <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_Af_fixed, test="Chisq")
test_fixed_4_5_uncertainty <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_C_fixed, test="Chisq")
test_fixed_4_5_agency <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_Ac_fixed, test="Chisq")
# covariates
test_fixed_4_5_pc <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_pc_fixed, test="Chisq")
test_fixed_4_5_wob <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_wob_fixed, test="Chisq")
test_fixed_4_5_gender <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_gender_fixed, test="Chisq")

## 6-7
test_fixed_6_7_EV <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_Af_fixed, test="Chisq")
test_fixed_6_7_uncertainty <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_C_fixed, test="Chisq")
test_fixed_6_7_agency <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_Ac_fixed, test="Chisq")
# covariates
test_fixed_6_7_pc <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_pc_fixed, test="Chisq")
test_fixed_6_7_wob <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_wob_fixed, test="Chisq")
test_fixed_6_7_gender <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_gender_fixed, test="Chisq")

## 8-9
test_fixed_8_9_EV <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_Af_fixed, test="Chisq")
test_fixed_8_9_uncertainty <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_C_fixed, test="Chisq")
test_fixed_8_9_agency <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_Ac_fixed, test="Chisq")
# covariates
test_fixed_8_9_pc <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_pc_fixed, test="Chisq")
test_fixed_8_9_wob <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_wob_fixed, test="Chisq")
test_fixed_8_9_gender <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_gender_fixed, test="Chisq")

## 10-12
test_fixed_10_12_EV <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_Af_fixed, test="Chisq")
test_fixed_10_12_uncertainty <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_C_fixed, test="Chisq")
test_fixed_10_12_agency <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_Ac_fixed, test="Chisq")
# covariates
test_fixed_10_12_pc <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_pc_fixed, test="Chisq")
test_fixed_10_12_wob <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_wob_fixed, test="Chisq")
test_fixed_10_12_gender <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_gender_fixed, test="Chisq")

#### Betas plots #### 
# Create betas plots
# Divided by age group
compt_hedonic_beta_by_group <- c(summary(compt_mod_4_5_full)$coefficients[2], summary(compt_mod_6_7_full)$coefficients[2], summary(compt_mod_8_9_full)$coefficients[2], summary(compt_mod_10_12_full)$coefficients[2], summary(adu_compt_mod_full)$coefficients[2])
compt_cognitive_beta_by_group <- c(summary(compt_mod_4_5_full)$coefficients[3], summary(compt_mod_6_7_full)$coefficients[3], summary(compt_mod_8_9_full)$coefficients[3], summary(compt_mod_10_12_full)$coefficients[3], summary(adu_compt_mod_full)$coefficients[3])
compt_instrumental_beta_by_group <- c(summary(compt_mod_4_5_full)$coefficients[4], summary(compt_mod_6_7_full)$coefficients[4], summary(compt_mod_8_9_full)$coefficients[4], summary(compt_mod_10_12_full)$coefficients[4], summary(adu_compt_mod_full)$coefficients[4])
compt_intercept_beta_by_group <- c(summary(compt_mod_4_5_full)$coefficients[1], summary(compt_mod_6_7_full)$coefficients[1], summary(compt_mod_8_9_full)$coefficients[1], summary(compt_mod_10_12_full)$coefficients[1], summary(adu_compt_mod_full)$coefficients[1])
compt_hedonic_sem_lower_by_group <- c(summary(compt_mod_4_5_full)$coefficients[2] - summary(compt_mod_4_5_full)$coefficients[2,2],
                                      summary(compt_mod_6_7_full)$coefficients[2] - summary(compt_mod_6_7_full)$coefficients[2,2],
                                      summary(compt_mod_8_9_full)$coefficients[2] - summary(compt_mod_8_9_full)$coefficients[2,2],
                                      summary(compt_mod_10_12_full)$coefficients[2] - summary(compt_mod_10_12_full)$coefficients[2,2],
                                      summary(adu_compt_mod_full)$coefficients[2] - summary(adu_compt_mod_full)$coefficients[2,2])
compt_cognitive_sem_lower_by_group <- c(summary(compt_mod_4_5_full)$coefficients[3] - summary(compt_mod_4_5_full)$coefficients[3,2],
                                        summary(compt_mod_6_7_full)$coefficients[3] - summary(compt_mod_6_7_full)$coefficients[3,2],
                                        summary(compt_mod_8_9_full)$coefficients[3] - summary(compt_mod_8_9_full)$coefficients[3,2],
                                        summary(compt_mod_10_12_full)$coefficients[3] - summary(compt_mod_10_12_full)$coefficients[3,2],
                                        summary(adu_compt_mod_full)$coefficients[3] - summary(adu_compt_mod_full)$coefficients[3,2])
compt_instrumental_sem_lower_by_group <-c(summary(compt_mod_4_5_full)$coefficients[4] - summary(compt_mod_4_5_full)$coefficients[4,2],
                                          summary(compt_mod_6_7_full)$coefficients[4] - summary(compt_mod_6_7_full)$coefficients[4,2],
                                          summary(compt_mod_8_9_full)$coefficients[4] - summary(compt_mod_8_9_full)$coefficients[4,2],
                                          summary(compt_mod_10_12_full)$coefficients[4] - summary(compt_mod_10_12_full)$coefficients[4,2],
                                          summary(adu_compt_mod_full)$coefficients[4] - summary(adu_compt_mod_full)$coefficients[4,2])
compt_intercept_sem_lower_by_group <- c(summary(compt_mod_4_5_full)$coefficients[1] - summary(compt_mod_4_5_full)$coefficients[1,2],
                                        summary(compt_mod_6_7_full)$coefficients[1] - summary(compt_mod_6_7_full)$coefficients[1,2],
                                        summary(compt_mod_8_9_full)$coefficients[1] - summary(compt_mod_8_9_full)$coefficients[1,2],
                                        summary(compt_mod_10_12_full)$coefficients[1] - summary(compt_mod_10_12_full)$coefficients[1,2],
                                        summary(adu_compt_mod_full)$coefficients[1] - summary(adu_compt_mod_full)$coefficients[1,2])
compt_hedonic_sem_upper_by_group <- c(summary(compt_mod_4_5_full)$coefficients[2] + summary(compt_mod_4_5_full)$coefficients[2,2],
                                      summary(compt_mod_6_7_full)$coefficients[2] + summary(compt_mod_6_7_full)$coefficients[2,2],
                                      summary(compt_mod_8_9_full)$coefficients[2] + summary(compt_mod_8_9_full)$coefficients[2,2],
                                      summary(compt_mod_10_12_full)$coefficients[2] + summary(compt_mod_10_12_full)$coefficients[2,2],
                                      summary(adu_compt_mod_full)$coefficients[2] + summary(adu_compt_mod_full)$coefficients[2,2])
compt_cognitive_sem_upper_by_group <- c(summary(compt_mod_4_5_full)$coefficients[3] + summary(compt_mod_4_5_full)$coefficients[3,2],
                                        summary(compt_mod_6_7_full)$coefficients[3] + summary(compt_mod_6_7_full)$coefficients[3,2],
                                        summary(compt_mod_8_9_full)$coefficients[3] + summary(compt_mod_8_9_full)$coefficients[3,2],
                                        summary(compt_mod_10_12_full)$coefficients[3] + summary(compt_mod_10_12_full)$coefficients[3,2],
                                        summary(adu_compt_mod_full)$coefficients[3] + summary(adu_compt_mod_full)$coefficients[3,2])
compt_instrumental_sem_upper_by_group <- c(summary(compt_mod_4_5_full)$coefficients[4] + summary(compt_mod_4_5_full)$coefficients[4,2],
                                           summary(compt_mod_6_7_full)$coefficients[4] + summary(compt_mod_6_7_full)$coefficients[4,2],
                                           summary(compt_mod_8_9_full)$coefficients[4] + summary(compt_mod_8_9_full)$coefficients[4,2],
                                           summary(compt_mod_10_12_full)$coefficients[4] + summary(compt_mod_10_12_full)$coefficients[4,2],
                                           summary(adu_compt_mod_full)$coefficients[4] + summary(adu_compt_mod_full)$coefficients[4,2])
compt_intercept_sem_upper_by_group <- c(summary(compt_mod_4_5_full)$coefficients[1] + summary(compt_mod_4_5_full)$coefficients[1,2],
                                        summary(compt_mod_6_7_full)$coefficients[1] + summary(compt_mod_6_7_full)$coefficients[1,2],
                                        summary(compt_mod_8_9_full)$coefficients[1] + summary(compt_mod_8_9_full)$coefficients[1,2],
                                        summary(compt_mod_10_12_full)$coefficients[1] + summary(compt_mod_10_12_full)$coefficients[1,2],
                                        summary(adu_compt_mod_full)$coefficients[1] + summary(adu_compt_mod_full)$coefficients[1,2])
hedonic_color_by_group <- c("#E41A1C", "#E41A1C", "#E41A1C", "#E41A1C", "#E41A1C")
cognitive_color_by_group <- c("#377EB8", "#377EB8", "#377EB8", "#377EB8", "#377EB8")
instrumental_color_by_group <- c("#4DAF4A", "#4DAF4A", "#4DAF4A", "#4DAF4A", "#4DAF4A")
intercept_color_by_group <- c("#808080", "#808080", "#808080", "#808080", "#808080")
x_points <- c(0.1, 0.2, 0.3, 0.4, 0.5)
x_labels <- c("4-5", "6-7", "8-9", "10-12", "Adults")

compt_betas_dat_by_group <- data.frame(compt_hedonic_beta_by_group, compt_cognitive_beta_by_group, compt_instrumental_beta_by_group, compt_intercept_beta_by_group,
                                       compt_hedonic_sem_lower_by_group, compt_cognitive_sem_lower_by_group, compt_instrumental_sem_lower_by_group, compt_intercept_sem_lower_by_group,
                                       compt_hedonic_sem_upper_by_group, compt_cognitive_sem_upper_by_group, compt_instrumental_sem_upper_by_group, compt_intercept_sem_upper_by_group,
                                       hedonic_color_by_group, cognitive_color_by_group, instrumental_color_by_group, intercept_color_by_group,
                                       x_points, x_labels)
if (rescale_variables) {
  hedonic_limits <- c(-0.5, 2.25)
  hedonic_breaks <- seq(-0.5, 2.25, 0.5)
  cognitive_limits <- c(-0.5, 2.25)
  cognitive_breaks <- seq(-0.5, 2.25, 0.5)
  instrumental_limits <- c(-0.5, 2.25)
  instrumental_breaks <- seq(-0.5, 2.25, 0.5)
} else {
  hedonic_limits <- c(-0.0015, 0.09)
  hedonic_breaks <- seq(-0.0015, 0.09, 0.015)
  cognitive_limits <- c(-0.2, 0.8)
  cognitive_breaks <- seq(-0.2, 0.8, 0.2)
  instrumental_limits <- c(0, 0.03)
  instrumental_breaks <- seq(0, 0.03, 0.005)
}
compt_hedonic_betas_by_group <- ggplot(compt_betas_dat_by_group) + 
  geom_line(aes(y=compt_hedonic_beta_by_group, x=x_points, colour = hedonic_color_by_group), size=1)+
  geom_point(aes(y=compt_hedonic_beta_by_group, x=x_points, colour = hedonic_color_by_group), size = 4)+
  geom_ribbon(aes(ymin=compt_hedonic_sem_lower_by_group, ymax=compt_hedonic_sem_upper_by_group, x=x_points, fill = "band"), alpha = 0.3)+
  scale_colour_manual("",values=compt_betas_dat_by_group$hedonic_color_by_group)+
  scale_fill_manual("",values=compt_betas_dat_by_group$hedonic_color_by_group) +
  theme_classic()+
  theme(text = element_text(size=16), legend.position="none") +
  scale_x_continuous(labels=compt_betas_dat_by_group$x_labels) +
  scale_y_continuous(limits=hedonic_limits, breaks=hedonic_breaks) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  labs(x="Age group", y="Standardized beta EV") +
  geom_hline(yintercept=0)

compt_cognitive_betas_by_group <- ggplot(compt_betas_dat_by_group) + 
  geom_line(aes(y=compt_cognitive_beta_by_group, x=x_points, colour = cognitive_color_by_group), size=1)+
  geom_point(aes(y=compt_cognitive_beta_by_group, x=x_points, colour = cognitive_color_by_group), size = 4)+
  geom_ribbon(aes(ymin=compt_cognitive_sem_lower_by_group, ymax=compt_cognitive_sem_upper_by_group, x=x_points, fill = "band"), alpha = 0.3)+
  scale_colour_manual("",values=compt_betas_dat_by_group$cognitive_color_by_group)+
  scale_fill_manual("",values=compt_betas_dat_by_group$cognitive_color_by_group) +
  theme_classic()+
  theme(text = element_text(size=15), legend.position="none") +
  scale_x_continuous(labels=compt_betas_dat_by_group$x_labels) +
  scale_y_continuous(limits=cognitive_limits, breaks=cognitive_breaks) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  labs(x="Age group", y="Standardized beta Uncertainty") +
  geom_hline(yintercept=0)

compt_instrumental_betas_by_group <- ggplot(compt_betas_dat_by_group) + 
  geom_line(aes(y=compt_instrumental_beta_by_group, x=x_points, colour = instrumental_color_by_group), size=1)+
  geom_point(aes(y=compt_instrumental_beta_by_group, x=x_points, colour = instrumental_color_by_group), size = 4)+
  geom_ribbon(aes(ymin=compt_instrumental_sem_lower_by_group, ymax=compt_instrumental_sem_upper_by_group, x=x_points, fill = "band"), alpha = 0.3)+
  scale_colour_manual("",values=compt_betas_dat_by_group$instrumental_color_by_group)+
  scale_fill_manual("",values=compt_betas_dat_by_group$instrumental_color_by_group) +
  theme_classic()+
  theme(text = element_text(size=15), legend.position="none") +
  scale_x_continuous(labels=compt_betas_dat_by_group$x_labels) +
  scale_y_continuous(limits=instrumental_limits, breaks=instrumental_breaks) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  labs(x="Age group", y="Standardized beta Agency") +
  geom_hline(yintercept=0)

grid.arrange(compt_hedonic_betas_by_group, compt_cognitive_betas_by_group, compt_instrumental_betas_by_group, ncol=3)

#### Individual fits ####
exclude_outliers <- FALSE
outlier_sd <- 3 # after how many sds consider an outlier
outlier_iqr <- 1 # after how many iqrs consider an outlier
## Children
id <- c()
for (i in unique(dat_compt$subject_ID)) {id <- append(id, i)}
age_group <- c()
age <- c()
pc <- c()
wob <- c()
gender <- c()
coeff_intercept_sbj <- c()
coeff_EV_sbj<-c()
coeff_uncertainty_sbj <- c()
coeff_agency_sbj <- c()
AIC_sbj <- c()
BIC_sbj <- c()
exclude_ids <- c()

for (i in seq_along(unique(dat_compt$subject_ID))) {
  id_sbj <- unique(dat_compt$subject_ID)[i]
  # organize in a data frame for subject = id_sbj
  dat <- subset(dat_compt, dat_compt$subject_ID == id_sbj)
  
  # update model-independent values
  age_group[i] <- dat$age_group[1]
  age[i] <- dat$age_in_years_non_Z[1]
  pc[i] <- dat$percent_comprehension_non_Z[1]
  wob[i] <- dat$wob_non_Z[1]
  gender[i] <- dat$gender[1]
  
  # logistic regression
  dat$delta_EV <- scale(dat$delta_EV_non_Z, center=FALSE)
  dat$delta_uncertainty_level <- scale(dat$delta_uncertainty_level_non_Z, center=FALSE)
  dat$delta_agency <- scale(dat$delta_agency_non_Z, center=FALSE)
  
  model <- glm(data=dat, info_choice ~ delta_EV + 
                 delta_uncertainty_level + 
                 delta_agency, family = "binomial", singular.ok=TRUE)
  
  coeff_intercept_sbj[i] <- model$coefficients[1]
  coeff_EV_sbj[i] <- model$coefficients["delta_EV"]
  coeff_uncertainty_sbj[i] <- model$coefficients["delta_uncertainty_level"]
  coeff_agency_sbj[i] <- model$coefficients["delta_agency"]
  
  if (!model$converged) {
    # print(paste(i, id_sbj, dat$age_group[1]))
    # print(paste("EV coeff:", round(coeff_EV_sbj[i], 2), 
    #             "uncertainty coeff:", round(coeff_uncertainty_sbj[i], 2),
    #             "agency coeff:", round(coeff_agency_sbj[i], 2)))
    # print(paste("EV ok:", ! glm(data=dat, info_choice ~ 
    #                               delta_uncertainty_level + 
    #                               delta_agency, family = "binomial", singular.ok=TRUE)$converged))
    # print(paste("uncertainty ok:", ! glm(data=dat, info_choice ~ delta_EV + 
    #                                        delta_agency, family = "binomial", singular.ok=TRUE)$converged))
    # print(paste("agency ok:", ! glm(data=dat, info_choice ~ delta_EV + 
    #                                   delta_uncertainty_level, 
    #                                 family = "binomial", singular.ok=TRUE)$converged))
    exclude_ids <- append(exclude_ids, id_sbj)
  }
  
  
  AIC_sbj[i]<- AIC(model)
  BIC_sbj[i]<- BIC(model)
  
}

# create dataframe
dat_compt_coeffs <- data.frame(id, age_group, age, pc, wob, gender, 
                               coeff_intercept_sbj, coeff_EV_sbj, 
                               coeff_uncertainty_sbj, coeff_agency_sbj, 
                               AIC_sbj, BIC_sbj)
dat_compt_coeffs$age_group <- factor(dat_compt_coeffs$age_group, 
                                     labels=c("4-5", "6-7", "8-9", "10-12"))
dat_compt_coeffs$age_group_coded <- "None"
dat_compt_coeffs$age_group_coded[dat_compt_coeffs$age  %in% c(4, 5)] <- 0
dat_compt_coeffs$age_group_coded[dat_compt_coeffs$age  %in% c(6, 7)] <- 1
dat_compt_coeffs$age_group_coded[dat_compt_coeffs$age  %in% c(8, 9)] <- 2
dat_compt_coeffs$age_group_coded[dat_compt_coeffs$age  %in% c(10, 11, 12)] <- 3
#dat_compt_coeffs$age_group_coded <- factor(dat_compt_coeffs$age_group_coded, labels=c("4-5", "6-7", "8-9","10-12"))
table(dat_compt_coeffs$age_group_coded)
#contrasts(dat_compt_coeffs$age_group_coded) <- contr.helmert(4)
dat_compt_coeffs$gender_coded <- as.factor(ifelse(dat_compt_coeffs$gender == "Male", 1, -1))
contrasts(dat_compt_coeffs$gender_coded) <- contr.helmert(2)
# exclude participants for which the model did not converge
dat_compt_coeffs <- subset(dat_compt_coeffs, ! id %in% exclude_ids)

## Adults
id <- c()
for (i in unique(adu_dat_compt$subject_ID)) {id <- append(id, i)}
age_group <- c()
age <- c()
pc <- c()
wob <- c()
gender <- c()
coeff_intercept_sbj <- c()
coeff_EV_sbj<-c()
coeff_uncertainty_sbj <- c()
coeff_agency_sbj <- c()
AIC_sbj <- c()
BIC_sbj <- c()
exclude_ids <- c()

for (i in seq_along(unique(adu_dat_compt$subject_ID))) {
  id_sbj <- unique(adu_dat_compt$subject_ID)[i]
  # organize in a data frame for subject = id_sbj
  dat <- subset(adu_dat_compt, adu_dat_compt$subject_ID == id_sbj)
  
  # update model-independent values
  age_group[i] <- dat$age_group[1]
  age[i] <- dat$age_in_years_non_Z[1]
  pc[i] <- dat$percent_comprehension_non_Z[1]
  wob[i] <- dat$wob_non_Z[1]
  gender[i] <- dat$gender[1]
  
  # logistic regression
  dat$delta_EV <- scale(dat$delta_EV_non_Z, center=FALSE)
  dat$delta_uncertainty_level <- scale(dat$delta_uncertainty_level_non_Z, center=FALSE)
  dat$delta_agency <- scale(dat$delta_agency_non_Z, center=FALSE)
  
  model <- glm(data=dat, info_choice ~ delta_EV + 
                 delta_uncertainty_level + 
                 delta_agency, family = "binomial", singular.ok=TRUE)
  
  coeff_intercept_sbj[i] <- model$coefficients[1]
  coeff_EV_sbj[i] <- model$coefficients["delta_EV"]
  coeff_uncertainty_sbj[i] <- model$coefficients["delta_uncertainty_level"]
  coeff_agency_sbj[i] <- model$coefficients["delta_agency"]
  
  if (!model$converged) {
    exclude_ids <- append(exclude_ids, id_sbj)
  }
  
  
  AIC_sbj[i]<- AIC(model)
  BIC_sbj[i]<- BIC(model)
  
}

# create dataframe
adu_dat_compt_coeffs <- data.frame(id, age_group, age, pc, wob, gender, 
                                   coeff_intercept_sbj, coeff_EV_sbj, 
                                   coeff_uncertainty_sbj, coeff_agency_sbj, 
                                   AIC_sbj, BIC_sbj)
adu_dat_compt_coeffs$age_group_coded <- 4
#contrasts(adu_dat_compt_coeffs$age_group_coded) <- contr.helmert(4)
adu_dat_compt_coeffs$gender_coded <- as.factor(ifelse(adu_dat_compt_coeffs$gender == "Male", 1, -1))
table(adu_dat_compt_coeffs$gender_coded)
contrasts(adu_dat_compt_coeffs$gender_coded) <- contr.helmert(2)
# exclude participants for which the model did not converge
adu_dat_compt_coeffs <- subset(adu_dat_compt_coeffs, ! id %in% exclude_ids)

# Merged dataframe
cols <- c("id", "age_group", "age", "pc", "wob", "gender", 
          "coeff_intercept_sbj", "coeff_EV_sbj", 
          "coeff_uncertainty_sbj", "coeff_agency_sbj", 
          "AIC_sbj", "BIC_sbj", "age_group_coded")
dat_compt_coeffs_child_adu <- rbind(dplyr::select(dat_compt_coeffs, all_of(cols)), 
                                    dplyr::select(adu_dat_compt_coeffs, all_of(cols)))


## Plots
compt_coeffs_EV_subj_violin_plot_by_group <- ggplot(dat_compt_coeffs_child_adu, aes(age_group_coded, coeff_EV_sbj)) +
  geom_violin(fill="#E41A1C", alpha=0.2) + 
  geom_point(color="#E41A1C") + 
  geom_hline(yintercept=0, color="black") +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.25, color = "black") +
  theme_classic() + 
  labs(x="Age group", y="Delta EV coefficient") 

# plot uncertainty coefficients
compt_coeffs_uncertainty_subj_violin_plot_by_group <- ggplot(dat_compt_coeffs_child_adu, aes(age_group_coded, coeff_uncertainty_sbj)) +
  geom_violin(fill="#377EB8", alpha=0.2) + 
  geom_point(color="#377EB8") + 
  geom_hline(yintercept=0, color="black") +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.25, color = "black") +
  theme_classic() + 
  labs(x="Age group", y="Delta uncertainty coefficient") 

# plot agency coefficients
compt_coeffs_agency_subj_violin_plot_by_group <- ggplot(dat_compt_coeffs_child_adu, aes(age_group_coded, coeff_agency_sbj)) +
  geom_violin(fill="#4DAF4A", alpha=0.2) + 
  geom_point(color="#4DAF4A") + 
  geom_hline(yintercept=0, color="black") +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.25, color = "black") +
  theme_classic() + 
  labs(x="Age group", y="Delta agency coefficient") 


grid.arrange(compt_coeffs_EV_subj_violin_plot_by_group, 
             compt_coeffs_uncertainty_subj_violin_plot_by_group, 
             compt_coeffs_agency_subj_violin_plot_by_group, ncol=3)

# T-tests 
test_EV_individual_4_5 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 0)$coeff_EV_sbj, mu=0)
test_uncertainty_individual_4_5 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 0)$coeff_uncertainty_sbj, mu=0)
test_agency_individual_4_5 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 0)$coeff_agency_sbj, mu=0)

test_EV_individual_6_7 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 1)$coeff_EV_sbj, mu=0)
test_uncertainty_individual_6_7 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 1)$coeff_uncertainty_sbj, mu=0)
test_agency_individual_6_7 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 2)$coeff_agency_sbj, mu=0)

test_EV_individual_8_9 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 2)$coeff_EV_sbj, mu=0)
test_uncertainty_individual_8_9 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 2)$coeff_uncertainty_sbj, mu=0)
test_agency_individual_8_9 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 2)$coeff_agency_sbj, mu=0)

test_EV_individual_10_12 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 3)$coeff_EV_sbj, mu=0)
test_uncertainty_individual_10_12 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 3)$coeff_uncertainty_sbj, mu=0)
test_agency_individual_10_12 <- t.test(subset(dat_compt_coeffs, dat_compt_coeffs$age_group_coded == 3)$coeff_agency_sbj, mu=0)

adu_test_EV_individual <- t.test(adu_dat_compt_coeffs$coeff_EV_sbj, mu=0)
adu_test_uncertainty_individual <- t.test(adu_dat_compt_coeffs$coeff_uncertainty_sbj, mu=0)
adu_test_agency_individual <- t.test(adu_dat_compt_coeffs$coeff_agency_sbj, mu=0)

