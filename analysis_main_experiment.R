# AUTHOR: GAIA MOLINARO
# Data analysis for the main experiment
# in Molinaro, Cogliati Dezza, & Sharot (in prep.)

# SET WORKING DIRECTORY
wd = "insert_your_wd_here"

#### LOAD REQUIRED PACKAGES #####
packages <- c("ggplot2", "GGally", "ggpubr", "afex", "plyr", "psych", "nls2", 
              "grid", "simr", "optimx", "plotrix", "BayesFactor", "ggcorrplot", 
              "tidyverse", "glue", "sjPlot", "sjmisc", "gridExtra", "rsq", 
              "sets", "qpcR")
pacman::p_load(packages, character.only = TRUE)

#### MAIN EXPERIMENT (COMPETITION) ####
## Load children's data
dat_compt <- read.csv("data_main_experiment_children.csv")
dat_compt <- subset(dat_compt, !(condition %in% c("catch_1", "catch_2")))
dat_compt <- subset(dat_compt, !(age_in_years == 13))
dat_compt <- subset(dat_compt, catch_trials_score == 100)
dat_compt <- subset(dat_compt, wob > 0.5)

# Create age groups
dat_compt$age_group <- "None"
dat_compt$age_group[dat_compt$age_in_years  %in% c(4, 5)] <- "4-5"
dat_compt$age_group[dat_compt$age_in_years  %in% c(6, 7)] <- "6-7"
dat_compt$age_group[dat_compt$age_in_years  %in% c(8, 9)] <- "8-9"
dat_compt$age_group[dat_compt$age_in_years  %in% c(10, 11, 12)] <- "10-12"
dat_compt$age_group <- factor(dat_compt$age_group, levels=c("4-5", "6-7", "8-9","10-12"))

# Turn variables into appropriate data types
# Use scale() for standardized values
dat_compt$delta_EV_non_Z <- as.numeric(dat_compt$delta_EV)
dat_compt$delta_uncertainty_level_non_Z <- as.numeric(dat_compt$delta_uncertainty_level)
dat_compt$delta_agency_non_Z <- as.numeric(dat_compt$delta_agency)
dat_compt$delta_SD <- as.numeric(dat_compt$delta_SD)
dat_compt$age_in_years_non_Z <- as.numeric(dat_compt$age_in_years)
dat_compt$age_in_years <- scale(as.numeric(dat_compt$age_in_years))
dat_compt$age_group_coded <- as.factor(dat_compt$age_group_coded)
dat_compt$RT_info_choice <- scale(as.numeric(dat_compt$RT_info_choice))
dat_compt$subject_ID <- as.factor(dat_compt$subject_ID)
dat_compt$gender_coded[dat_compt$gender_coded==3] <- 1 # 1 = male, 2 = female, 3 = other
dat_compt$gender_coded <- as.factor(dat_compt$gender_coded)
contrasts(dat_compt$gender_coded) <- contr.helmert(2)
dat_compt$info_choice <- as.factor(dat_compt$info_choice) # 0 = left, 1 = right
dat_compt$percent_comprehension_non_Z <- as.numeric(dat_compt$percent_comprehension)
dat_compt$percent_comprehension <- scale(as.numeric(dat_compt$percent_comprehension))
dat_compt$wob_non_Z <- as.numeric(dat_compt$wob)
dat_compt$wob <- scale(as.numeric(dat_compt$wob))
dat_compt$chance <- 0.5
# correct fishing choices based on non-decoy items
dat_compt$fishing_w_info_non_Z_1 <- as.numeric(dat_compt$prop_correct_fishing_with_info_1)
dat_compt$fishing_w_info_1 <- scale(as.numeric(dat_compt$prop_correct_fishing_with_info_1))
dat_compt$fishing_no_info_non_Z_1 <- as.numeric(dat_compt$prop_correct_fishing_without_info_1)
dat_compt$fishing_no_info_1 <- scale(as.numeric(dat_compt$prop_correct_fishing_without_info_1))
dat_compt$fishing_non_Z_1 <- as.numeric(dat_compt$prop_correct_fishing_1)
dat_compt$fishing_1 <- scale(as.numeric(dat_compt$prop_correct_fishing_1))
# correct fishing choices based on what participants see
dat_compt$fishing_w_info_non_Z_2 <- as.numeric(dat_compt$prop_correct_fishing_with_info_2)
dat_compt$fishing_w_info_2 <- scale(as.numeric(dat_compt$prop_correct_fishing_with_info_2))
dat_compt$fishing_no_info_non_Z_2 <- as.numeric(dat_compt$prop_correct_fishing_without_info_2)
dat_compt$fishing_no_info_2 <- scale(as.numeric(dat_compt$prop_correct_fishing_without_info_2))
dat_compt$fishing_non_Z_2 <- as.numeric(dat_compt$prop_correct_fishing_2)
dat_compt$fishing_2 <- scale(as.numeric(dat_compt$prop_correct_fishing_2))


## Compute scaled deltas
# Find min and max and rescale between -1 and 1
rescale_deltas <- function(dat, deltas=c("EV", "uncertainty_level", "agency")){
  if ("EV" %in% deltas){
    min_EV = min(as.numeric(c(dat$EV_L, dat$EV_R)))
    max_EV = max(as.numeric(c(dat$EV_L, dat$EV_R)))
    dat$EV_L <- 2*((as.numeric(dat$EV_L) - min_EV)/(max_EV - min_EV)) -1
    dat$EV_R <- 2*((as.numeric(dat$EV_R) - min_EV)/(max_EV - min_EV)) -1
    dat$delta_EV <- dat$EV_R - dat$EV_L
    }
  if ("uncertainty_level" %in% deltas) {
    min_uncertainty_level = min(as.numeric(c(dat$uncertainty_level_L, dat$uncertainty_level_R)))
  max_uncertainty_level = max(as.numeric(c(dat$uncertainty_level_L, dat$uncertainty_level_R)))
  dat$uncertainty_level_L <- 2*((as.numeric(dat$uncertainty_level_L) - min_uncertainty_level)/(max_uncertainty_level - min_uncertainty_level)) -1
  dat$uncertainty_level_R <- 2*((as.numeric(dat$uncertainty_level_R) - min_uncertainty_level)/(max_uncertainty_level - min_uncertainty_level)) -1
  dat$delta_uncertainty_level <- dat$uncertainty_level_R - dat$uncertainty_level_L
  }
  if ("agency" %in% deltas) {
    min_agency = min(as.numeric(c(dat$agency_probL, dat$agency_probR)))
  max_agency = max(as.numeric(c(dat$agency_probL, dat$agency_probR)))
  dat$agency_probL <- 2*((as.numeric(dat$agency_probL) - min_agency)/(max_agency - min_agency)) -1
  dat$agency_probR <- 2*((as.numeric(dat$agency_probR) - min_agency)/(max_agency - min_agency)) -1
  dat$delta_agency <- dat$agency_probR - dat$agency_probL
  }
  return(dat)
}

dat_compt <- rescale_deltas(dat_compt)

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
dat_compt_4_5$percent_comprehension <- scale(dat_compt_4_5$percent_comprehension_non_Z)
dat_compt_6_7$percent_comprehension <- scale(dat_compt_6_7$percent_comprehension_non_Z)
dat_compt_8_9$percent_comprehension <- scale(dat_compt_8_9$percent_comprehension_non_Z)
dat_compt_10_12$percent_comprehension <- scale(dat_compt_10_12$percent_comprehension_non_Z)
dat_compt_4_5$wob <- scale(dat_compt_4_5$wob_non_Z)
dat_compt_6_7$wob <- scale(dat_compt_6_7$wob_non_Z)
dat_compt_8_9$wob <- scale(dat_compt_8_9$wob_non_Z)
dat_compt_10_12$wob <- scale(dat_compt_10_12$wob_non_Z)
dat_compt_4_5 <- rescale_deltas(dat_compt_4_5)
dat_compt_6_7 <- rescale_deltas(dat_compt_6_7)
dat_compt_8_9 <- rescale_deltas(dat_compt_8_9)
dat_compt_10_12 <- rescale_deltas(dat_compt_10_12)

## Load adults' data
setwd(wd)
adu_dat_compt <- read.csv("data_main_experiment_adults.csv")
adu_dat_compt <- subset(adu_dat_compt, !(condition %in% c("catch_1", "catch_2")))
adu_dat_compt <- subset(adu_dat_compt, (catch_trials_score == 100) & (wob > 0.5))

# Fix trial number
adu_dat_compt$trial_number <- substr(adu_dat_compt$condition, 7, 8)

# Turn variables into appropriate data types
adu_dat_compt$delta_EV_non_Z <- as.numeric(adu_dat_compt$delta_EV)
adu_dat_compt$delta_uncertainty_level_non_Z <- as.numeric(adu_dat_compt$delta_uncertainty_level)
adu_dat_compt$delta_agency_non_Z <- as.numeric(adu_dat_compt$delta_agency)
adu_dat_compt$delta_SD <- as.numeric(adu_dat_compt$delta_SD)
adu_dat_compt$age_in_years_non_Z <- as.numeric(adu_dat_compt$age_in_years)
adu_dat_compt$age_in_years <- scale(as.numeric(adu_dat_compt$age_in_years))
adu_dat_compt$age_group <- "Adult"
adu_dat_compt$age_group_coded <- as.factor(4)
adu_dat_compt$RT_info_choice <- scale(as.numeric(adu_dat_compt$RT_info_choice))
adu_dat_compt$subject_ID <- as.factor(adu_dat_compt$subject_ID)
adu_dat_compt$gender_coded[adu_dat_compt$gender_coded==3] <- 1 # contrast will be female vs non-female
adu_dat_compt$gender_coded <- as.factor(adu_dat_compt$gender_coded)
contrasts(adu_dat_compt$gender_coded) <- contr.helmert(2)
adu_dat_compt$info_choice <- as.factor(adu_dat_compt$info_choice) # 0 = left, 1 = right
adu_dat_compt$percent_comprehension_non_Z <- as.numeric(adu_dat_compt$percent_comprehension)
adu_dat_compt$percent_comprehension <- scale(as.numeric(adu_dat_compt$percent_comprehension))
adu_dat_compt$wob_non_Z <- as.numeric(adu_dat_compt$wob)
adu_dat_compt$wob <- scale(as.numeric(adu_dat_compt$wob))
adu_dat_compt$chance <- 0.5
# correct fishing choices based on non-decoy items
adu_dat_compt$fishing_w_info_non_Z_1 <- as.numeric(adu_dat_compt$prop_correct_fishing_with_info_1)
adu_dat_compt$fishing_w_info_1 <- scale(as.numeric(adu_dat_compt$prop_correct_fishing_with_info_1))
adu_dat_compt$fishing_no_info_non_Z_1 <- as.numeric(adu_dat_compt$prop_correct_fishing_without_info_1)
adu_dat_compt$fishing_no_info_1 <- scale(as.numeric(adu_dat_compt$prop_correct_fishing_without_info_1))
adu_dat_compt$fishing_non_Z_1 <- as.numeric(adu_dat_compt$prop_correct_fishing_1)
adu_dat_compt$fishing_1 <- scale(as.numeric(adu_dat_compt$prop_correct_fishing_1))
# correct fishing choices based on what participants see
adu_dat_compt$fishing_w_info_non_Z_2 <- as.numeric(adu_dat_compt$prop_correct_fishing_with_info_2)
adu_dat_compt$fishing_w_info_2 <- scale(as.numeric(adu_dat_compt$prop_correct_fishing_with_info_2))
adu_dat_compt$fishing_no_info_non_Z_2 <- as.numeric(adu_dat_compt$prop_correct_fishing_without_info_2)
adu_dat_compt$fishing_no_info_2 <- scale(as.numeric(adu_dat_compt$prop_correct_fishing_without_info_2))
adu_dat_compt$fishing_non_Z_2 <- as.numeric(adu_dat_compt$prop_correct_fishing_2)
adu_dat_compt$fishing_2 <- scale(as.numeric(adu_dat_compt$prop_correct_fishing_2))


## Compute scaled deltas
# Find min and max and rescale between -1 and 1
adu_dat_compt <- rescale_deltas(adu_dat_compt)

# Check rescaled and raw deltas have the same sign (adults)
check_delta_signs(adu_dat_compt)

## Merged datasets (children + adults)
# Create data subsets
cols <- c("subject_ID", "info_choice", 
          "delta_EV", "delta_uncertainty_level", "delta_agency", 
          "EV_L", "EV_R", "uncertainty_level_L", "uncertainty_level_R",
          "agency_probL", "agency_probR",
          "gender_coded", "wob", "percent_comprehension", 
          "wob_non_Z", "percent_comprehension_non_Z", 
          "age_group", "age_group_coded",
          "fishing_w_info_non_Z_1", "fishing_no_info_non_Z_1", "fishing_non_Z_1",
          "fishing_w_info_non_Z_2", "fishing_no_info_non_Z_2", "fishing_non_Z_2")
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


# Re-Z-score variables
dat_compt_child_adu$percent_comprehension <- scale(dat_compt_child_adu$percent_comprehension_non_Z)
dat_compt_4_5_adu$percent_comprehension <- scale(dat_compt_4_5_adu$percent_comprehension_non_Z)
dat_compt_6_7_adu$percent_comprehension <- scale(dat_compt_6_7_adu$percent_comprehension_non_Z)
dat_compt_8_9_adu$percent_comprehension <- scale(dat_compt_8_9_adu$percent_comprehension_non_Z)
dat_compt_10_12_adu$percent_comprehension <- scale(dat_compt_10_12_adu$percent_comprehension_non_Z)
dat_compt_child_adu$wob <- scale(dat_compt_child_adu$wob_non_Z)
dat_compt_4_5_adu$wob <- scale(dat_compt_4_5_adu$wob_non_Z)
dat_compt_6_7_adu$wob <- scale(dat_compt_6_7_adu$wob_non_Z)
dat_compt_8_9_adu$wob <- scale(dat_compt_8_9_adu$wob_non_Z)
dat_compt_10_12_adu$wob <- scale(dat_compt_10_12_adu$wob_non_Z)
dat_compt_child_adu <- rescale_deltas(dat_compt_child_adu)
dat_compt_4_5_adu <- rescale_deltas(dat_compt_4_5_adu)
dat_compt_6_7_adu <- rescale_deltas(dat_compt_6_7_adu)
dat_compt_8_9_adu <- rescale_deltas(dat_compt_8_9_adu)
dat_compt_10_12_adu <- rescale_deltas(dat_compt_10_12_adu)

# Set gender contrasts
contrasts(dat_compt_child_adu$gender_coded) <- contr.helmert(2)
contrasts(dat_compt_4_5_adu$gender_coded) <- contr.helmert(2)
contrasts(dat_compt_6_7_adu$gender_coded) <- contr.helmert(2)
contrasts(dat_compt_8_9_adu$gender_coded) <- contr.helmert(2)
contrasts(dat_compt_10_12_adu$gender_coded) <- contr.helmert(2)

# Create new variable and  contrasts for children vs adults factor
# somehow
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

# Within group
t.test(subset(dat_compt, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(adu_dat_compt, trial_number == 1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_compt_4_5, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_compt_6_7, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_compt_8_9, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_compt_10_12, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)


# Between groups
t.test(subset(dat_compt, trial_number ==1)$percent_comprehension_non_Z, 
       subset(adu_dat_compt, trial_number == 1)$percent_comprehension_non_Z, var.equal=TRUE)
t.test(subset(dat_compt, trial_number ==1)$wob_non_Z, 
       subset(adu_dat_compt, trial_number == 1)$wob_non_Z, var.equal=TRUE)

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
# Children vs adults
compt_pc <- ggplot(subset(dat_compt_child_adu_pc_wob, variable2 %in% 
  c("percent_comprehension_child", "percent_comprehension_adu")), 
  aes(x=variable2, y=score, fill=variable2)) +
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
compt_pc


compt_wob <- ggplot(subset(dat_compt_child_adu_pc_wob, variable2 %in% 
                            c("wob_child", "wob_adu")), 
                   aes(x=variable2, y=score, fill=variable2)) +
  geom_jitter(height=0, color = subset(dat_compt_child_adu_pc_wob, 
                                       variable2 %in% c("wob_child", "wob_adu"))$color,
              size=1.5, alpha=0.9)+   
  geom_violin(alpha=0.3) +
  geom_hline(yintercept=0.5, color="#AAAAAA") +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0, 1)) + 
  scale_fill_manual(values=c("#AAAAAA", "#444444")) +
  #stat_summary(fun=mean, geom="point", shape=20, size=4, color="black", fill="black") +
  theme_classic() + labs(x="comprehension", y="Score") + theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm")) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))
compt_wob

grid.arrange(compt_pc, compt_wob, ncol=2)

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
compt_pc_by_group

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
compt_wob_by_group

grid.arrange(compt_pc_by_group, compt_wob_by_group, ncol=2)

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
# Children overall 
# Chance only
child_compt_mod_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID),
                                data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, with interactions
child_compt_mod_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + gender_coded + age_in_years +
                                age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                (delta_EV + delta_uncertainty_level + delta_agency
                                 | subject_ID), 
                              data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, with previous trial as a predictor

child_compt_mod_full_pers <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + gender_coded + age_in_years + prev_info_choice +
                                age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                (delta_EV + delta_uncertainty_level + delta_agency
                                 | subject_ID), 
                              data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)



# Full, without intercept
child_compt_mod_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + gender_coded + age_in_years - 1 +
                                age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                (delta_EV + delta_uncertainty_level + delta_agency
                                 | subject_ID), 
                              data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full,  without interaction between delta_agency and age
child_compt_mod_full_drop_Af_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                            percent_comprehension + wob + gender_coded + age_in_years +
                                            age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                            (delta_EV + delta_uncertainty_level + delta_agency
                                             | subject_ID), 
                                          data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without interaction between delta_uncertainty_level and age
child_compt_mod_full_drop_C_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                           percent_comprehension + wob + gender_coded + age_in_years +
                                           age_in_years:delta_EV + age_in_years:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency
                                            | subject_ID), 
                                         data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without interaction between delta_agency and age
child_compt_mod_full_drop_Ac_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                            percent_comprehension + wob + gender_coded + age_in_years +
                                            age_in_years:delta_EV + age_in_years:delta_uncertainty_level +
                                            (delta_EV + delta_uncertainty_level + delta_agency
                                             | subject_ID), 
                                          data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without EV as a fixed effect
child_compt_mod_full_drop_Af_fixed <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + gender_coded + age_in_years +
                                age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                (delta_EV + delta_uncertainty_level + delta_agency
                                 | subject_ID), 
                              data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without uncertainty as a fixed effect
child_compt_mod_full_drop_C_fixed <- glmer(info_choice ~ delta_EV + delta_agency +
                                              percent_comprehension + wob + gender_coded + age_in_years +
                                              age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without agency as a fixed effect
child_compt_mod_full_drop_Ac_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                              percent_comprehension + wob + gender_coded + age_in_years +
                                              age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without comprehension as a fixed effect
child_compt_mod_full_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               wob + gender_coded + age_in_years +
                                              age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without wob as a fixed effect
child_compt_mod_full_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              percent_comprehension + gender_coded + age_in_years +
                                              age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without gender as a fixed effect
child_compt_mod_full_drop_gender_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob  + age_in_years +
                                              age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, without age as a fixed effect
child_compt_mod_full_drop_age_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob + gender_coded +
                                              age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, without EV as a fixed or interaction effect
child_compt_mod_full_drop_Af_fixed_and_int <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob + gender_coded + age_in_years +
                                              age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without EV as a fixed or interaction effect
child_compt_mod_full_drop_C_fixed_and_int <- glmer(info_choice ~ delta_EV + delta_agency +
                                                      percent_comprehension + wob + gender_coded + age_in_years +
                                                      age_in_years:delta_EV + age_in_years:delta_agency +
                                                      (delta_EV + delta_uncertainty_level + delta_agency
                                                       | subject_ID), 
                                                    data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without agency as a fixed or interaction effect
child_compt_mod_full_drop_Ac_fixed_and_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                                      percent_comprehension + wob + gender_coded + age_in_years +
                                                      age_in_years:delta_EV + age_in_years:delta_uncertainty_level +
                                                      (delta_EV + delta_uncertainty_level + delta_agency
                                                       | subject_ID), 
                                                    data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Affect only
child_compt_mod_Af <- glmer(info_choice ~ delta_EV +
                              percent_comprehension + wob + gender_coded + age_in_years +
                              delta_EV:age_in_years +
                              (delta_EV 
                               | subject_ID), 
                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition only 
child_compt_mod_C <- glmer(info_choice ~  delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + age_in_years +
                             delta_uncertainty_level:age_in_years +
                             (delta_uncertainty_level 
                              | subject_ID),
                           data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Action only
child_compt_mod_Ac <- glmer(info_choice ~ delta_agency +
                              percent_comprehension + wob + gender_coded + age_in_years +
                              delta_agency:age_in_years +
                              (delta_agency 
                               | subject_ID), 
                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Cognition 
child_compt_mod_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               percent_comprehension + wob + gender_coded + age_in_years +
                               delta_EV:age_in_years + delta_uncertainty_level:age_in_years +
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID),
                             data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Action 
child_compt_mod_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                percent_comprehension + wob + gender_coded + age_in_years +
                                delta_EV:age_in_years + delta_agency:age_in_years +
                                (delta_EV + delta_agency
                                 | subject_ID),
                              data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition +  Action
child_compt_mod_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + age_in_years +
                               delta_uncertainty_level:age_in_years + delta_agency:age_in_years +
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID),
                             data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

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
# Full, with previous trial as a predictor
compt_mod_4_5_full_pers <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                              percent_comprehension + wob + gender_coded + prev_info_choice +
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

# 6-7
# Chance only
compt_mod_6_7_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                              data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
                             
# Full 
compt_mod_6_7_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_uncertainty_level + delta_agency                               | subject_ID), 
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

#### Significance tests within groups ####
# (within single groups)
# access results of test a with a$Chisq, a$Df, a$Pr[2]
## Children overall
# delta EV: overall, only fixed, fixed + interaction
anova(child_compt_mod_full, child_compt_mod_CAc, test="Chisq")
child_test_fixed_EV <- anova(child_compt_mod_full, child_compt_mod_full_drop_Af_fixed, test="Chisq")
anova(child_compt_mod_full, child_compt_mod_full_drop_Af_fixed_and_int, test="Chisq")
# delta uncertainty: overall, only fixed, fixed + interaction
anova(child_compt_mod_full, child_compt_mod_AfAc, test="Chisq")
child_test_fixed_uncertainty <- anova(child_compt_mod_full, child_compt_mod_full_drop_C_fixed, test="Chisq")
anova(child_compt_mod_full, child_compt_mod_full_drop_C_fixed_and_int, test="Chisq")
# delta agency: overall, only fixed, fixed + interaction
anova(child_compt_mod_full, child_compt_mod_AfC, test="Chisq")
child_test_fixed_agency <-anova(child_compt_mod_full, child_compt_mod_full_drop_Ac_fixed, test="Chisq")
anova(child_compt_mod_full, child_compt_mod_full_drop_Ac_fixed_and_int, test="Chisq")
# delta EV x age
child_test_fixed_EV_int <-anova(child_compt_mod_full, child_compt_mod_full_drop_Af_int, test="Chisq")
# delta uncertainty x age
child_test_fixed_uncertainty_int <- anova(child_compt_mod_full, child_compt_mod_full_drop_C_int, test="Chisq")
# delta agency x age
child_test_fixed_agency_int <- anova(child_compt_mod_full, child_compt_mod_full_drop_Ac_int, test="Chisq")
# covariates
child_test_fixed_pc <- anova(child_compt_mod_full, child_compt_mod_full_drop_pc_fixed, test="Chisq")
child_test_fixed_wob <- anova(child_compt_mod_full, child_compt_mod_full_drop_wob_fixed, test="Chisq")
child_test_fixed_gender <- anova(child_compt_mod_full, child_compt_mod_full_drop_gender_fixed, test="Chisq")
child_test_fixed_age <- anova(child_compt_mod_full, child_compt_mod_full_drop_age_fixed, test="Chisq")

## Adults overall
# delta EV: overall, only fixed
anova(adu_compt_mod_full, adu_compt_mod_CAc, test="Chisq")
adu_test_fixed_EV <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_Af_fixed, test="Chisq")
# delta uncertainty: overall, only fixed
anova(adu_compt_mod_full, adu_compt_mod_AfAc, test="Chisq")
adu_test_fixed_uncertainty <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_C_fixed, test="Chisq")
# delta agency: overall, only fixed
anova(adu_compt_mod_full, adu_compt_mod_AfC, test="Chisq")
adu_test_fixed_agency <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_Ac_fixed, test="Chisq")
# covariates
adu_test_fixed_pc <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_pc_fixed, test="Chisq")
adu_test_fixed_wob <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_wob_fixed, test="Chisq")
adu_test_fixed_gender <- anova(adu_compt_mod_full, adu_compt_mod_full_drop_gender_fixed, test="Chisq")

## 4-5
# delta EV: overall, only fixed
anova(compt_mod_4_5_full, compt_mod_4_5_CAc, test="Chisq")
test_fixed_4_5_EV <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_Af_fixed, test="Chisq")
# delta uncertainty: overall, only fixed
anova(compt_mod_4_5_full, compt_mod_4_5_AfAc, test="Chisq")
test_fixed_4_5_uncertainty <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_C_fixed, test="Chisq")
# delta agency: overall, only fixed
anova(compt_mod_4_5_full, compt_mod_4_5_AfC, test="Chisq")
test_fixed_4_5_agency <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_Ac_fixed, test="Chisq")
# covariates
test_fixed_4_5_pc <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_pc_fixed, test="Chisq")
test_fixed_4_5_wob <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_wob_fixed, test="Chisq")
test_fixed_4_5_gender <- anova(compt_mod_4_5_full, compt_mod_4_5_full_drop_gender_fixed, test="Chisq")
# additional test: previous trial
test_fixed_4_5_prev_trial <- anova(compt_mod_4_5_full_pers, compt_mod_4_5_full, test="Chisq")


## 6-7
# delta EV: overall, only fixed
anova(compt_mod_6_7_full, compt_mod_6_7_CAc, test="Chisq")
test_fixed_6_7_EV <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_Af_fixed, test="Chisq")
# delta uncertainty: overall, only fixed
anova(compt_mod_6_7_full, compt_mod_6_7_AfAc, test="Chisq")
test_fixed_6_7_uncertainty <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_C_fixed, test="Chisq")
# delta agency: overall, only fixed
anova(compt_mod_6_7_full, compt_mod_6_7_AfC, test="Chisq")
test_fixed_6_7_agency <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_Ac_fixed, test="Chisq")
# covariates
test_fixed_6_7_pc <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_pc_fixed, test="Chisq")
test_fixed_6_7_wob <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_wob_fixed, test="Chisq")
test_fixed_6_7_gender <- anova(compt_mod_6_7_full, compt_mod_6_7_full_drop_gender_fixed, test="Chisq")


## 8-9
# delta EV: overall, only fixed
anova(compt_mod_8_9_full, compt_mod_8_9_CAc, test="Chisq")
test_fixed_8_9_EV <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_Af_fixed, test="Chisq")
# delta uncertainty: overall, only fixed
anova(compt_mod_8_9_full, compt_mod_8_9_AfAc, test="Chisq")
test_fixed_8_9_uncertainty <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_C_fixed, test="Chisq")
# delta agency: overall, only fixed
anova(compt_mod_8_9_full, compt_mod_8_9_AfC, test="Chisq")
test_fixed_8_9_agency <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_Ac_fixed, test="Chisq")
# covariates
test_fixed_8_9_pc <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_pc_fixed, test="Chisq")
test_fixed_8_9_wob <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_wob_fixed, test="Chisq")
test_fixed_8_9_gender <- anova(compt_mod_8_9_full, compt_mod_8_9_full_drop_gender_fixed, test="Chisq")


## 10-12
# delta EV: overall, only fixed
anova(compt_mod_10_12_full, compt_mod_10_12_CAc, test="Chisq")
test_fixed_10_12_EV <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_Af_fixed, test="Chisq")
# delta uncertainty: overall, only fixed
anova(compt_mod_10_12_full, compt_mod_10_12_AfAc, test="Chisq")
test_fixed_10_12_uncertainty <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_C_fixed, test="Chisq")
# delta agency: overall, only fixed
anova(compt_mod_10_12_full, compt_mod_10_12_AfC, test="Chisq")
test_fixed_10_12_agency <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_Ac_fixed, test="Chisq")
# covariates
test_fixed_10_12_pc <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_pc_fixed, test="Chisq")
test_fixed_10_12_wob <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_wob_fixed, test="Chisq")
test_fixed_10_12_gender <- anova(compt_mod_10_12_full, compt_mod_10_12_full_drop_gender_fixed, test="Chisq")


# Create csv file with tests of fixed effects
variable <- c("children_EV",  "children_uncertainty", "children_agency", 
              "children_comprehension", "children_wob", "children_gender", 
              "children_age", "children_age_by_EV", 
              "children_age_by_uncertainty", "children_age_by_agency",
              "adults_EV",  "adults_uncertainty", "adults_agency", 
              "adults_comprehension", "adults_wob", "adults_gender", 
              "4_5_EV",  "4_5_uncertainty", "4_5_agency", 
              "4_5_comprehension", "4_5_wob", "4_5_gender", 
              "6_7_EV",  "6_7_uncertainty", "6_7_agency", 
              "6_7_comprehension", "6_7_wob", "6_7_gender", 
              "8_9_EV",  "8_9_uncertainty", "8_9_agency", 
              "8_9_comprehension", "8_9_wob", "8_9_gender", 
              "10_12_EV",  "10_12_uncertainty", "10_12_agency", 
              "10_12_comprehension", "10_12_wob", "10_12_gender")

# store coefficients
child_compt_mod_full_coeffs <- summary(child_compt_mod_full)$coefficients[, 1]
adu_compt_mod_full_coeffs <- summary(adu_compt_mod_full)$coefficients[, 1]
compt_mod_4_5_full_coeffs <- summary(compt_mod_4_5_full)$coefficients[, 1]
compt_mod_6_7_full_coeffs <- summary(compt_mod_6_7_full)$coefficients[, 1]
compt_mod_8_9_full_coeffs <- summary(compt_mod_8_9_full)$coefficients[, 1]
compt_mod_10_12_full_coeffs <- summary(compt_mod_10_12_full)$coefficients[, 1]

# store sems
child_compt_mod_full_sems <- summary(child_compt_mod_full)$coefficients[, 2]
adu_compt_mod_full_sems <- summary(adu_compt_mod_full)$coefficients[, 2]
compt_mod_4_5_full_sems <- summary(compt_mod_4_5_full)$coefficients[, 2]
compt_mod_6_7_full_sems <- summary(compt_mod_6_7_full)$coefficients[, 2]
compt_mod_8_9_full_sems <- summary(compt_mod_8_9_full)$coefficients[, 2]
compt_mod_10_12_full_sems <- summary(compt_mod_10_12_full)$coefficients[, 2]


estimate <- c(child_compt_mod_full_coeffs["delta_EV"], 
              child_compt_mod_full_coeffs["delta_uncertainty_level"],
              child_compt_mod_full_coeffs["delta_agency"],
              child_compt_mod_full_coeffs["percent_comprehension"],
              child_compt_mod_full_coeffs["wob"],
              child_compt_mod_full_coeffs["gender_coded1"],
              child_compt_mod_full_coeffs["age_in_years"],
              child_compt_mod_full_coeffs["delta_EV:age_in_years"],
              child_compt_mod_full_coeffs["delta_uncertainty_level:age_in_years"],
              child_compt_mod_full_coeffs["delta_agency:age_in_years"],
              adu_compt_mod_full_coeffs["delta_EV"], 
              adu_compt_mod_full_coeffs["delta_uncertainty_level"],
              adu_compt_mod_full_coeffs["delta_agency"],
              adu_compt_mod_full_coeffs["percent_comprehension"],
              adu_compt_mod_full_coeffs["wob"],
              adu_compt_mod_full_coeffs["gender_coded1"],
              compt_mod_4_5_full_coeffs["delta_EV"], 
              compt_mod_4_5_full_coeffs["delta_uncertainty_level"],
              compt_mod_4_5_full_coeffs["delta_agency"],
              compt_mod_4_5_full_coeffs["percent_comprehension"],
              compt_mod_4_5_full_coeffs["wob"],
              compt_mod_4_5_full_coeffs["gender_coded1"],
              compt_mod_6_7_full_coeffs["delta_EV"], 
              compt_mod_6_7_full_coeffs["delta_uncertainty_level"],
              compt_mod_6_7_full_coeffs["delta_agency"],
              compt_mod_6_7_full_coeffs["percent_comprehension"],
              compt_mod_6_7_full_coeffs["wob"],
              compt_mod_6_7_full_coeffs["gender_coded1"],
              compt_mod_8_9_full_coeffs["delta_EV"], 
              compt_mod_8_9_full_coeffs["delta_uncertainty_level"],
              compt_mod_8_9_full_coeffs["delta_agency"],
              compt_mod_8_9_full_coeffs["percent_comprehension"],
              compt_mod_8_9_full_coeffs["wob"],
              compt_mod_8_9_full_coeffs["gender_coded1"],
              compt_mod_10_12_full_coeffs["delta_EV"], 
              compt_mod_10_12_full_coeffs["delta_uncertainty_level"],
              compt_mod_10_12_full_coeffs["delta_agency"],
              compt_mod_10_12_full_coeffs["percent_comprehension"],
              compt_mod_10_12_full_coeffs["wob"],
              compt_mod_10_12_full_coeffs["gender_coded1"])

sem <- c(child_compt_mod_full_sems["delta_EV"], 
              child_compt_mod_full_sems["delta_uncertainty_level"],
              child_compt_mod_full_sems["delta_agency"],
              child_compt_mod_full_sems["percent_comprehension"],
              child_compt_mod_full_sems["wob"],
              child_compt_mod_full_sems["gender_coded1"],
              child_compt_mod_full_sems["age_in_years"],
              child_compt_mod_full_sems["delta_EV:age_in_years"],
              child_compt_mod_full_sems["delta_uncertainty_level:age_in_years"],
              child_compt_mod_full_sems["delta_agency:age_in_years"],
              adu_compt_mod_full_sems["delta_EV"], 
              adu_compt_mod_full_sems["delta_uncertainty_level"],
              adu_compt_mod_full_sems["delta_agency"],
              adu_compt_mod_full_sems["percent_comprehension"],
              adu_compt_mod_full_sems["wob"],
              adu_compt_mod_full_sems["gender_coded1"],
              compt_mod_4_5_full_sems["delta_EV"], 
              compt_mod_4_5_full_sems["delta_uncertainty_level"],
              compt_mod_4_5_full_sems["delta_agency"],
              compt_mod_4_5_full_sems["percent_comprehension"],
              compt_mod_4_5_full_sems["wob"],
              compt_mod_4_5_full_sems["gender_coded1"],
              compt_mod_6_7_full_sems["delta_EV"], 
              compt_mod_6_7_full_sems["delta_uncertainty_level"],
              compt_mod_6_7_full_sems["delta_agency"],
              compt_mod_6_7_full_sems["percent_comprehension"],
              compt_mod_6_7_full_sems["wob"],
              compt_mod_6_7_full_sems["gender_coded1"],
              compt_mod_8_9_full_sems["delta_EV"], 
              compt_mod_8_9_full_sems["delta_uncertainty_level"],
              compt_mod_8_9_full_sems["delta_agency"],
              compt_mod_8_9_full_sems["percent_comprehension"],
              compt_mod_8_9_full_sems["wob"],
              compt_mod_8_9_full_sems["gender_coded1"],
              compt_mod_10_12_full_sems["delta_EV"], 
              compt_mod_10_12_full_sems["delta_uncertainty_level"],
              compt_mod_10_12_full_sems["delta_agency"],
              compt_mod_10_12_full_sems["percent_comprehension"],
              compt_mod_10_12_full_sems["wob"],
              compt_mod_10_12_full_sems["gender_coded1"])


chisq <- c(child_test_fixed_EV$Chisq[2], child_test_fixed_uncertainty$Chisq[2],
           child_test_fixed_agency$Chisq[2], child_test_fixed_pc$Chisq[2],
           child_test_fixed_wob$Chisq[2], child_test_fixed_gender$Chisq[2], 
           child_test_fixed_age$Chisq[2], child_test_fixed_EV_int$Chisq[2], 
           child_test_fixed_uncertainty_int$Chisq[2], child_test_fixed_agency_int$Chisq[2],
           adu_test_fixed_EV$Chisq[2], adu_test_fixed_uncertainty$Chisq[2],
           adu_test_fixed_agency$Chisq[2], adu_test_fixed_pc$Chisq[2],
           adu_test_fixed_wob$Chisq[2], adu_test_fixed_gender$Chisq[2],
           test_fixed_4_5_EV$Chisq[2], test_fixed_4_5_uncertainty$Chisq[2],
           test_fixed_4_5_agency$Chisq[2], test_fixed_4_5_pc$Chisq[2],
           test_fixed_4_5_wob$Chisq[2], test_fixed_4_5_gender$Chisq[2],
           test_fixed_6_7_EV$Chisq[2], test_fixed_6_7_uncertainty$Chisq[2],
           test_fixed_6_7_agency$Chisq[2], test_fixed_6_7_pc$Chisq[2],
           test_fixed_6_7_wob$Chisq[2], test_fixed_6_7_gender$Chisq[2],
           test_fixed_8_9_EV$Chisq[2], test_fixed_8_9_uncertainty$Chisq[2],
           test_fixed_8_9_agency$Chisq[2], test_fixed_8_9_pc$Chisq[2],
           test_fixed_8_9_wob$Chisq[2], test_fixed_8_9_gender$Chisq[2],
           test_fixed_10_12_EV$Chisq[2], test_fixed_10_12_uncertainty$Chisq[2],
           test_fixed_10_12_agency$Chisq[2], test_fixed_10_12_pc$Chisq[2],
           test_fixed_10_12_wob$Chisq[2], test_fixed_10_12_gender$Chisq[2])

df <- c(child_test_fixed_EV$Df[2], child_test_fixed_uncertainty$Df[2],
           child_test_fixed_agency$Df[2], child_test_fixed_pc$Df[2],
           child_test_fixed_wob$Df[2], child_test_fixed_gender$Df[2], 
           child_test_fixed_age$Df[2], child_test_fixed_EV_int$Df[2], 
           child_test_fixed_uncertainty_int$Df[2], child_test_fixed_agency_int$Df[2],
           adu_test_fixed_EV$Df[2], adu_test_fixed_uncertainty$Df[2],
           adu_test_fixed_agency$Df[2], adu_test_fixed_pc$Df[2],
           adu_test_fixed_wob$Df[2], adu_test_fixed_gender$Df[2],
           test_fixed_4_5_EV$Df[2], test_fixed_4_5_uncertainty$Df[2],
           test_fixed_4_5_agency$Df[2], test_fixed_4_5_pc$Df[2],
           test_fixed_4_5_wob$Df[2], test_fixed_4_5_gender$Df[2],
           test_fixed_6_7_EV$Df[2], test_fixed_6_7_uncertainty$Df[2],
           test_fixed_6_7_agency$Df[2], test_fixed_6_7_pc$Df[2],
           test_fixed_6_7_wob$Df[2], test_fixed_6_7_gender$Df[2],
           test_fixed_8_9_EV$Df[2], test_fixed_8_9_uncertainty$Df[2],
           test_fixed_8_9_agency$Df[2], test_fixed_8_9_pc$Df[2],
           test_fixed_8_9_wob$Df[2], test_fixed_8_9_gender$Df[2],
           test_fixed_10_12_EV$Df[2], test_fixed_10_12_uncertainty$Df[2],
           test_fixed_10_12_agency$Df[2], test_fixed_10_12_pc$Df[2],
           test_fixed_10_12_wob$Df[2], test_fixed_10_12_gender$Df[2])

p_value <- c(child_test_fixed_EV$Pr[2], child_test_fixed_uncertainty$Pr[2],
        child_test_fixed_agency$Pr[2], child_test_fixed_pc$Pr[2],
        child_test_fixed_wob$Pr[2], child_test_fixed_gender$Pr[2], 
        child_test_fixed_age$Pr[2], child_test_fixed_EV_int$Pr[2], 
        child_test_fixed_uncertainty_int$Pr[2], child_test_fixed_agency_int$Pr[2],
        adu_test_fixed_EV$Pr[2], adu_test_fixed_uncertainty$Pr[2],
        adu_test_fixed_agency$Pr[2], adu_test_fixed_pc$Pr[2],
        adu_test_fixed_wob$Pr[2], adu_test_fixed_gender$Pr[2],
        test_fixed_4_5_EV$Pr[2], test_fixed_4_5_uncertainty$Pr[2],
        test_fixed_4_5_agency$Pr[2], test_fixed_4_5_pc$Pr[2],
        test_fixed_4_5_wob$Pr[2], test_fixed_4_5_gender$Pr[2],
        test_fixed_6_7_EV$Pr[2], test_fixed_6_7_uncertainty$Pr[2],
        test_fixed_6_7_agency$Pr[2], test_fixed_6_7_pc$Pr[2],
        test_fixed_6_7_wob$Pr[2], test_fixed_6_7_gender$Pr[2],
        test_fixed_8_9_EV$Pr[2], test_fixed_8_9_uncertainty$Pr[2],
        test_fixed_8_9_agency$Pr[2], test_fixed_8_9_pc$Pr[2],
        test_fixed_8_9_wob$Pr[2], test_fixed_8_9_gender$Pr[2],
        test_fixed_10_12_EV$Pr[2], test_fixed_10_12_uncertainty$Pr[2],
        test_fixed_10_12_agency$Pr[2], test_fixed_10_12_pc$Pr[2],
        test_fixed_10_12_wob$Pr[2], test_fixed_10_12_gender$Pr[2])


compt_mods_fixed_effects <- data.frame(variable, estimate, sem, chisq, df, p_value)

# Save fixed effects as a file
write.csv(compt_mods_fixed_effects, sprintf("%s/compt_mods_fixed_effects.csv", wd), row.names = FALSE)

#### AIC #### 
# Create data frame with AIC scores
age_group <- factor(c("4-5", "6-7", "8-9", "10-12", "Adults"), levels=c("4-5", "6-7", "8-9", "10-12", "Adults"))
chance_AIC <- c(AIC(compt_mod_4_5_chance), AIC(compt_mod_6_7_chance), AIC(compt_mod_8_9_chance), AIC(compt_mod_10_12_chance), AIC(adu_compt_mod_chance))
Af_AIC <- c(AIC(compt_mod_4_5_Af), AIC(compt_mod_6_7_Af), AIC(compt_mod_8_9_Af), AIC(compt_mod_10_12_Af), AIC(adu_compt_mod_Af))
C_AIC <- c(AIC(compt_mod_4_5_C), AIC(compt_mod_6_7_C), AIC(compt_mod_8_9_C), AIC(compt_mod_10_12_C), AIC(adu_compt_mod_C))
Ac_AIC <- c(AIC(compt_mod_4_5_Ac), AIC(compt_mod_6_7_Ac), AIC(compt_mod_8_9_Ac), AIC(compt_mod_10_12_Ac), AIC(adu_compt_mod_Ac))
AfC_AIC <- c(AIC(compt_mod_4_5_AfC), AIC(compt_mod_6_7_AfC), AIC(compt_mod_8_9_AfC), AIC(compt_mod_10_12_AfC), AIC(adu_compt_mod_AfC))
AfAc_AIC <- c(AIC(compt_mod_4_5_AfAc), AIC(compt_mod_6_7_AfAc), AIC(compt_mod_8_9_AfAc), AIC(compt_mod_10_12_AfAc), AIC(adu_compt_mod_AfAc))
CAc_AIC <- c(AIC(compt_mod_4_5_CAc), AIC(compt_mod_6_7_CAc),AIC(compt_mod_8_9_CAc),AIC(compt_mod_10_12_CAc), AIC(adu_compt_mod_CAc))
full_AIC <- c(AIC(compt_mod_4_5_full), AIC(compt_mod_6_7_full), AIC(compt_mod_8_9_full), AIC(compt_mod_10_12_full), AIC(adu_compt_mod_full))

AIC_scores <- data.frame(age_group, chance_AIC,
                         Af_AIC, C_AIC, Ac_AIC, AfC_AIC, AfAc_AIC, CAc_AIC, full_AIC)

# With covariates and random effects
write.csv(AIC_scores,sprintf("%s/AIC_scores.csv", wd), row.names = FALSE)

AIC_scores_long <- gather(AIC_scores, model, AIC_score, chance_AIC:full_AIC, factor_key=TRUE)
AIC_scores_long$model <- factor(AIC_scores_long$model, labels=c("Chance", 
                         "EV", "Uncertainty", "Agency",
                        "EV + Uncertainty", "EV + Agency", "Uncertainty + Agency", "EV + Uncertainty + Agency"))
AIC_scores_plot <- ggplot(data=AIC_scores_long, aes(x=model, y=AIC_score)) +
  geom_bar(stat="identity") + 
  facet_grid(~age_group) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=60, hjust=1)) +
  theme(axis.text.y = element_text(size=14))
AIC_scores_plot

# Children vs adults
age_group <- factor(c("Children", "Adults"), levels=c("Children", "Adults"))
chance_AIC <- c(AIC(child_compt_mod_chance), AIC(adu_compt_mod_chance))
Af_AIC <- c(AIC(child_compt_mod_Af), AIC(adu_compt_mod_Af))
C_AIC <- c(AIC(child_compt_mod_C), AIC(adu_compt_mod_C))
Ac_AIC <- c(AIC(child_compt_mod_Ac), AIC(adu_compt_mod_Ac))
AfC_AIC <- c(AIC(child_compt_mod_AfC), AIC(adu_compt_mod_AfC))
AfAc_AIC <- c(AIC(child_compt_mod_AfAc), AIC(adu_compt_mod_AfAc))
CAc_AIC <- c(AIC(child_compt_mod_CAc), AIC(adu_compt_mod_CAc))
full_AIC <- c(AIC(child_compt_mod_full), AIC(adu_compt_mod_full))

AIC_scores_compt_child_vs_adu <- data.frame(age_group, chance_AIC, 
                                            Af_AIC, C_AIC, Ac_AIC, AfC_AIC, AfAc_AIC, CAc_AIC, full_AIC)
# With covariates and random effects
write.csv(AIC_scores_compt_child_vs_adu, sprintf("%s/AIC_scores_compt_child_vs_adu.csv", wd), row.names = FALSE)

AIC_scores_compt_child_vs_adu_long <- gather(AIC_scores_compt_child_vs_adu, model, AIC_score, chance_AIC:full_AIC, factor_key=TRUE)
AIC_scores_compt_child_vs_adu_long$model <- factor(AIC_scores_compt_child_vs_adu_long$model, labels=c("Chance", 
                                                                                                      "EV", "Uncertainty", "Agency",
                                                                                                            "EV + Uncertainty", "EV + Agency", "Uncertainty + Agency", "EV + Uncertainty + Agency"))
AIC_scores_compt_child_vs_adu_plot <- ggplot(data=AIC_scores_compt_child_vs_adu_long, aes(x=model, y=AIC_score)) +
  geom_bar(stat="identity") + 
  facet_grid(~age_group) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=60, hjust=1, size=14)) +
  theme(axis.text.y = element_text(size=14))
AIC_scores_compt_child_vs_adu_plot

#### R^2 #### 
# Create data frame with R^2 scores (full models only)
age_group <- factor(c("4-5", "6-7", "8-9", "10-12", "Adults"), levels=c("4-5", "6-7", "8-9", "10-12", "Adults"))
x_points <- c(0.1, 0.2, 0.3, 0.4, 0.5)
x_labels <- c("4-5", "6-7", "8-9", "10-12", "Adults")
R2 <- c(rsq(compt_mod_4_5_full)$model, rsq(compt_mod_6_7_full)$model, rsq(compt_mod_8_9_full)$model, rsq(compt_mod_10_12_full)$model, rsq(adu_compt_mod_full)$model)
one_min_R2 <- 1-R2
R2_scores <- data.frame(age_group, R2, one_min_R2, x_points, x_labels)
R2_pol_1 = data.frame(x_R2_pol_1=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.4, 0.3, 0.2, 0.1), y_R2_pol_1=c(0, 0, 0, 0, 0, one_min_R2[5],one_min_R2[4], one_min_R2[3], one_min_R2[2], one_min_R2[1]))
R2_pol_2 = data.frame(x_R2_pol_2=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.4, 0.3, 0.2, 0.1), y_R2_pol_2=c(one_min_R2[1],one_min_R2[2], one_min_R2[3], one_min_R2[4], one_min_R2[5], 1, 1, 1, 1, 1))
R2_pol_3 = data.frame(x_R2_pol_3=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.4, 0.3, 0.2, 0.1), y_R2_pol_3=c(0, 0, 0, 0, 0, R2[5],R2[4], R2[3], R2[2], R2[1]))
R2_pol_4 = data.frame(x_R2_pol_4=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.4, 0.3, 0.2, 0.1), y_R2_pol_4=c(R2[1],R2[2], R2[3], R2[4], R2[5], 1, 1, 1, 1, 1))


one_min_R2_scores_plot <- ggplot(R2_scores) + 
  geom_polygon(data=R2_pol_1, aes(x=x_R2_pol_1, y=y_R2_pol_1), fill = "dark grey") +
  geom_polygon(data=R2_pol_2, aes(x=x_R2_pol_2, y=y_R2_pol_2), fill = "light grey") +
  geom_line(aes(y=one_min_R2, x=x_points), size=1) +
  geom_point(aes(y=one_min_R2, x=x_points), size = 4) +
  theme_classic() + theme(text = element_text(size=15), legend.position="none") +
  scale_x_continuous(labels=R2_scores$x_labels) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  #geom_hline(yintercept=1) +
  #geom_vline(xintercept=c(0.1, 0.5)) +
  labs(x="Age group", y="Noise\n(1- R^2)") 
one_min_R2_scores_plot 

R2_scores_plot <- ggplot(R2_scores) + 
  #geom_polygon(data=R2_pol_3, aes(x=x_R2_pol_3, y=y_R2_pol_3), fill = "light grey") +
  #geom_polygon(data=R2_pol_4, aes(x=x_R2_pol_4, y=y_R2_pol_4), fill = "dark grey") +
  geom_line(aes(y=R2, x=x_points), size=1) +
  geom_point(aes(y=R2, x=x_points), size = 4) +
  theme_classic() + theme(text = element_text(size=15), legend.position="none") +
  scale_x_continuous(labels=R2_scores$x_labels) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  #geom_hline(yintercept=1) +
  #geom_vline(xintercept=c(0.1, 0.5)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  labs(x="Age group", y="R^2") 
R2_scores_plot 

#### Betas plots #### 
# Create betas plots
# Children vs adults
compt_hedonic_beta <- c(summary(child_compt_mod_full)$coefficients[2], summary(adu_compt_mod_full)$coefficients[2])
compt_cognitive_beta <- c(summary(child_compt_mod_full)$coefficients[3], summary(adu_compt_mod_full)$coefficients[3])
compt_instrumental_beta <- c(summary(child_compt_mod_full)$coefficients[4], summary(adu_compt_mod_full)$coefficients[4])
compt_hedonic_sem_lower <- c(summary(child_compt_mod_full)$coefficients[2] - summary(child_compt_mod_full)$coefficients[2,2],
                             summary(adu_compt_mod_full)$coefficients[2] - summary(adu_compt_mod_full)$coefficients[2,2])
compt_cognitive_sem_lower <- c(summary(child_compt_mod_full)$coefficients[3] - summary(child_compt_mod_full)$coefficients[3,2],
                             summary(adu_compt_mod_full)$coefficients[3] - summary(adu_compt_mod_full)$coefficients[3,2])
compt_instrumental_sem_lower <- c(summary(child_compt_mod_full)$coefficients[4] - summary(child_compt_mod_full)$coefficients[4,2],
                               summary(adu_compt_mod_full)$coefficients[4] - summary(adu_compt_mod_full)$coefficients[4,2])
compt_hedonic_sem_upper <- c(summary(child_compt_mod_full)$coefficients[2] + summary(child_compt_mod_full)$coefficients[2,2],
                             summary(adu_compt_mod_full)$coefficients[2] + summary(adu_compt_mod_full)$coefficients[2,2])
compt_cognitive_sem_upper <- c(summary(child_compt_mod_full)$coefficients[3] + summary(child_compt_mod_full)$coefficients[3,2],
                               summary(adu_compt_mod_full)$coefficients[3] + summary(adu_compt_mod_full)$coefficients[3,2])
compt_instrumental_sem_upper <- c(summary(child_compt_mod_full)$coefficients[4] + summary(child_compt_mod_full)$coefficients[4,2],
                                  summary(adu_compt_mod_full)$coefficients[4] + summary(adu_compt_mod_full)$coefficients[4,2])
hedonic_color <- "#E41A1C"
cognitive_color <- "#377EB8"
instrumental_color <- "#4DAF4A"

compt_betas_dat <- data.frame(compt_hedonic_beta, compt_cognitive_beta, compt_instrumental_beta,
                              compt_hedonic_sem_lower, compt_cognitive_sem_lower, compt_instrumental_sem_lower,
                              compt_hedonic_sem_upper, compt_cognitive_sem_upper, compt_instrumental_sem_upper,
                              hedonic_color, cognitive_color, instrumental_color)

compt_motive_betas <- ggplot()+
  scale_y_continuous(breaks=seq(0,2,0.1), limits=c(0, 2), expand=c(0, 0)) + 
  # Children
  geom_errorbar(aes(0.3, ymin = compt_hedonic_sem_lower[1], ymax = compt_hedonic_sem_upper[1]), size = 0.8, width = .05, colour = hedonic_color) +
  geom_point(aes(0.3, compt_hedonic_beta[1]), size = 4, color = hedonic_color)+
  geom_errorbar(aes(0.6, ymin = compt_cognitive_sem_lower[1], ymax = compt_cognitive_sem_upper[1]), size = 0.8, width = .05, colour = cognitive_color) +
  geom_point(aes(0.6, compt_cognitive_beta[1]), size = 4, color = cognitive_color)+
  geom_errorbar(aes(0.9, ymin = compt_instrumental_sem_lower[1], ymax = compt_instrumental_sem_upper[1]), size = 0.8, width = .05, colour = instrumental_color) +
  geom_point(aes(0.9, compt_instrumental_beta[1]), size = 4, color = instrumental_color)+
  # Adults
  geom_errorbar(aes(0.4, ymin = compt_hedonic_sem_lower[2], ymax = compt_hedonic_sem_upper[2]), size = 0.8, width = .05, colour = hedonic_color, linetype="dashed") +
  geom_point(aes(0.4, compt_hedonic_beta[2]), size = 4, color = hedonic_color)+
  geom_errorbar(aes(0.7, ymin = compt_cognitive_sem_lower[2], ymax = compt_cognitive_sem_upper[2]), size = 0.8, width = .05, colour = cognitive_color, linetype="dashed") +
  geom_point(aes(0.7, compt_cognitive_beta[2]), size = 4, color = cognitive_color)+
  geom_errorbar(aes(1, ymin = compt_instrumental_sem_lower[2], ymax = compt_instrumental_sem_upper[2]), size = 0.8, width = .05, colour = instrumental_color, linetype="dashed") +
  geom_point(aes(1, compt_instrumental_beta[2]), size = 4, color = instrumental_color)+
  geom_vline(xintercept=c(0.5, 0.8), linetype="dotted") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size=20)) +
  labs(x="Utilities", y="Standardized beta predicting information-seeking")
compt_motive_betas

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

compt_hedonic_betas_by_group <- ggplot(compt_betas_dat_by_group) + 
  geom_line(aes(y=compt_hedonic_beta_by_group, x=x_points, colour = hedonic_color_by_group), size=1)+
  geom_point(aes(y=compt_hedonic_beta_by_group, x=x_points, colour = hedonic_color_by_group), size = 4)+
  geom_ribbon(aes(ymin=compt_hedonic_sem_lower_by_group, ymax=compt_hedonic_sem_upper_by_group, x=x_points, fill = "band"), alpha = 0.3)+
  scale_colour_manual("",values=compt_betas_dat_by_group$hedonic_color_by_group)+
  scale_fill_manual("",values=compt_betas_dat_by_group$hedonic_color_by_group) +
  theme_classic()+
  theme(text = element_text(size=16), legend.position="none") +
  scale_x_continuous(labels=compt_betas_dat_by_group$x_labels) +
  scale_y_continuous(limits=c(-0.5, 2.25), breaks=seq(-0.5, 2.25, 0.5)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  labs(x="Age group", y="Standardized beta EV") +
  geom_hline(yintercept=0)
  #+ geom_vline(xintercept=0.45, linetype="dotted")
compt_hedonic_betas_by_group

compt_cognitive_betas_by_group <- ggplot(compt_betas_dat_by_group) + 
  geom_line(aes(y=compt_cognitive_beta_by_group, x=x_points, colour = cognitive_color_by_group), size=1)+
  geom_point(aes(y=compt_cognitive_beta_by_group, x=x_points, colour = cognitive_color_by_group), size = 4)+
  geom_ribbon(aes(ymin=compt_cognitive_sem_lower_by_group, ymax=compt_cognitive_sem_upper_by_group, x=x_points, fill = "band"), alpha = 0.3)+
  scale_colour_manual("",values=compt_betas_dat_by_group$cognitive_color_by_group)+
  scale_fill_manual("",values=compt_betas_dat_by_group$cognitive_color_by_group) +
  theme_classic()+
  theme(text = element_text(size=15), legend.position="none") +
  scale_x_continuous(labels=compt_betas_dat_by_group$x_labels) +
  scale_y_continuous(limits=c(-0.5, 2.25), breaks=seq(-0.5, 2.25, 0.5)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  labs(x="Age group", y="Standardized beta Uncertainty") +
  geom_hline(yintercept=0)
  #+ geom_vline(xintercept=0.45, linetype="dotted")
compt_cognitive_betas_by_group

compt_instrumental_betas_by_group <- ggplot(compt_betas_dat_by_group) + 
  geom_line(aes(y=compt_instrumental_beta_by_group, x=x_points, colour = instrumental_color_by_group), size=1)+
  geom_point(aes(y=compt_instrumental_beta_by_group, x=x_points, colour = instrumental_color_by_group), size = 4)+
  geom_ribbon(aes(ymin=compt_instrumental_sem_lower_by_group, ymax=compt_instrumental_sem_upper_by_group, x=x_points, fill = "band"), alpha = 0.3)+
  scale_colour_manual("",values=compt_betas_dat_by_group$instrumental_color_by_group)+
  scale_fill_manual("",values=compt_betas_dat_by_group$instrumental_color_by_group) +
  theme_classic()+
  theme(text = element_text(size=15), legend.position="none") +
  scale_x_continuous(labels=compt_betas_dat_by_group$x_labels) +
  scale_y_continuous(limits=c(-0.5, 2.25), breaks=seq(-0.5, 2.25, 0.5)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  labs(x="Age group", y="Standardized beta Agency") +
  geom_hline(yintercept=0)
  #+ geom_vline(xintercept=0.45, linetype="dotted")
compt_instrumental_betas_by_group

compt_intercept_betas_by_group <- ggplot(compt_betas_dat_by_group) + 
  geom_line(aes(y=compt_intercept_beta_by_group, x=x_points, colour = intercept_color_by_group), size=1)+
  geom_point(aes(y=compt_intercept_beta_by_group, x=x_points, colour = intercept_color_by_group), size = 4)+
  geom_ribbon(aes(ymin=compt_intercept_sem_lower_by_group, ymax=compt_intercept_sem_upper_by_group, x=x_points, fill = "band"), alpha = 0.3)+
  scale_colour_manual("",values=compt_betas_dat_by_group$intercept_color_by_group)+
  scale_fill_manual("",values=compt_betas_dat_by_group$intercept_color_by_group) +
  theme_classic()+
  theme(text = element_text(size=15), legend.position="none") +
  scale_x_continuous(labels=compt_betas_dat_by_group$x_labels) +
  scale_y_continuous(limits=c(-1.75, 2.25), breaks=seq(-1.75, 2.25, 0.5)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  labs(x="Age group", y="Standardized beta coefficient Intercept") +
  geom_hline(yintercept=0)
#+ geom_vline(xintercept=0.45, linetype="dotted")
compt_intercept_betas_by_group

grid.arrange(compt_hedonic_betas_by_group, compt_cognitive_betas_by_group, compt_instrumental_betas_by_group, ncol=3)

####  Children vs adults tests #### 
## Models
# Children vs adults
# Full
compt_mod_child_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + group + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                           data = dat_compt_child_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and EV
compt_mod_child_adu_drop_Af_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + group + 
                               group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and uncertainty
compt_mod_child_adu_drop_C_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + group + 
                               group:delta_EV + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and agency
compt_mod_child_adu_drop_Ac_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + group + 
                               group:delta_EV + group:delta_uncertainty_level +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop EV fixed
compt_mod_child_adu_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + group + 
                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Drop uncertainty fixed
compt_mod_child_adu_drop_C_fixed <- glmer(info_choice ~ delta_EV + delta_agency +
                               percent_comprehension + wob + gender_coded + group + 
                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop agency fixed
compt_mod_child_adu_drop_Ac_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               percent_comprehension + wob + gender_coded + group + 
                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop comprehension fixed
compt_mod_child_adu_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                               wob + gender_coded + group + 
                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop wob fixed
compt_mod_child_adu_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                               percent_comprehension + gender_coded + group + 
                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop gender fixed 
compt_mod_child_adu_drop_gender_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + group + 
                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop group fixed
compt_mod_child_adu_drop_group_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + 
                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# 4-5 year-olds
# Full
compt_mod_4_5_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + group + 
                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_4_5_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and EV
compt_mod_4_5_adu_drop_Af_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                           percent_comprehension + wob + gender_coded + group + 
                                           group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_4_5_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and uncertainty
compt_mod_4_5_adu_drop_C_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                          percent_comprehension + wob + gender_coded + group + 
                                          group:delta_EV + group:delta_agency +
                                          (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                        data = dat_compt_4_5_adu, family = binomial, 
                                        control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and agency
compt_mod_4_5_adu_drop_Ac_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                           percent_comprehension + wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_4_5_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop EV fixed
compt_mod_4_5_adu_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                             percent_comprehension + wob + gender_coded + group + 
                                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                           data = dat_compt_4_5_adu, family = binomial, 
                                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop uncertainty fixed
compt_mod_4_5_adu_drop_C_fixed <- glmer(info_choice ~ delta_EV + delta_agency +
                                            percent_comprehension + wob + gender_coded + group + 
                                            group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                            (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                          data = dat_compt_4_5_adu, family = binomial, 
                                          control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop agency fixed
compt_mod_4_5_adu_drop_Ac_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                             percent_comprehension + wob + gender_coded + group + 
                                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                           data = dat_compt_4_5_adu, family = binomial, 
                                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop comprehension fixed
compt_mod_4_5_adu_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                             wob + gender_coded + group + 
                                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                           data = dat_compt_4_5_adu, family = binomial, 
                                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop wob fixed
compt_mod_4_5_adu_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              percent_comprehension + gender_coded + group + 
                                              group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                            data = dat_compt_4_5_adu, family = binomial, 
                                            control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop gender fixed 
compt_mod_4_5_adu_drop_gender_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                 percent_comprehension + wob + group + 
                                                 group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                               data = dat_compt_4_5_adu, family = binomial, 
                                               control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop group fixed
compt_mod_4_5_adu_drop_group_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + 
                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_compt_4_5_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# 6-7-year-olds
# Full
compt_mod_6_7_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + group + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                           data = dat_compt_6_7_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and EV
compt_mod_6_7_adu_drop_Af_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                         percent_comprehension + wob + gender_coded + group + 
                                         group:delta_uncertainty_level + group:delta_agency +
                                         (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                       data = dat_compt_6_7_adu, family = binomial, 
                                       control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and uncertainty
compt_mod_6_7_adu_drop_C_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                        percent_comprehension + wob + gender_coded + group + 
                                        group:delta_EV + group:delta_agency +
                                        (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                      data = dat_compt_6_7_adu, family = binomial, 
                                      control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and agency
compt_mod_6_7_adu_drop_Ac_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                         percent_comprehension + wob + gender_coded + group + 
                                         group:delta_EV + group:delta_uncertainty_level +
                                         (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                       data = dat_compt_6_7_adu, family = binomial, 
                                       control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop EV fixed
compt_mod_6_7_adu_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                           percent_comprehension + wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_6_7_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop uncertainty fixed
compt_mod_6_7_adu_drop_C_fixed <- glmer(info_choice ~ delta_EV + delta_agency +
                                          percent_comprehension + wob + gender_coded + group + 
                                          group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                          (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                        data = dat_compt_6_7_adu, family = binomial, 
                                        control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop agency fixed
compt_mod_6_7_adu_drop_Ac_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                           percent_comprehension + wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_6_7_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop comprehension fixed
compt_mod_6_7_adu_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                           wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_6_7_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop wob fixed
compt_mod_6_7_adu_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                            percent_comprehension + gender_coded + group + 
                                            group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                            (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                          data = dat_compt_6_7_adu, family = binomial, 
                                          control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop gender fixed 
compt_mod_6_7_adu_drop_gender_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob + group + 
                                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                             data = dat_compt_6_7_adu, family = binomial, 
                                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop group fixed
compt_mod_6_7_adu_drop_group_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                           data = dat_compt_6_7_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# 8-9-year-olds
# Full
compt_mod_8_9_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + group + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                           data = dat_compt_8_9_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and EV
compt_mod_8_9_adu_drop_Af_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                         percent_comprehension + wob + gender_coded + group + 
                                         group:delta_uncertainty_level + group:delta_agency +
                                         (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                       data = dat_compt_8_9_adu, family = binomial, 
                                       control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and uncertainty
compt_mod_8_9_adu_drop_C_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                        percent_comprehension + wob + gender_coded + group + 
                                        group:delta_EV + group:delta_agency +
                                        (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                      data = dat_compt_8_9_adu, family = binomial, 
                                      control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and agency
compt_mod_8_9_adu_drop_Ac_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                         percent_comprehension + wob + gender_coded + group + 
                                         group:delta_EV + group:delta_uncertainty_level +
                                         (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                       data = dat_compt_8_9_adu, family = binomial, 
                                       control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop EV fixed
compt_mod_8_9_adu_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                           percent_comprehension + wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_8_9_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop uncertainty fixed
compt_mod_8_9_adu_drop_C_fixed <- glmer(info_choice ~ delta_EV + delta_agency +
                                          percent_comprehension + wob + gender_coded + group + 
                                          group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                          (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                        data = dat_compt_8_9_adu, family = binomial, 
                                        control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop agency fixed
compt_mod_8_9_adu_drop_Ac_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                           percent_comprehension + wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_8_9_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop comprehension fixed
compt_mod_8_9_adu_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                           wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_8_9_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop wob fixed
compt_mod_8_9_adu_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                            percent_comprehension + gender_coded + group + 
                                            group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                            (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                          data = dat_compt_8_9_adu, family = binomial, 
                                          control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop gender fixed 
compt_mod_8_9_adu_drop_gender_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob + group + 
                                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                             data = dat_compt_8_9_adu, family = binomial, 
                                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop group fixed
compt_mod_8_9_adu_drop_group_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                           data = dat_compt_8_9_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# 10-12 year-olds
# Full
compt_mod_10_12_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + group + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                           data = dat_compt_10_12_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and EV
compt_mod_10_12_adu_drop_Af_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                         percent_comprehension + wob + gender_coded + group + 
                                         group:delta_uncertainty_level + group:delta_agency +
                                         (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                       data = dat_compt_10_12_adu, family = binomial, 
                                       control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and uncertainty
compt_mod_10_12_adu_drop_C_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                        percent_comprehension + wob + gender_coded + group + 
                                        group:delta_EV + group:delta_agency +
                                        (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                      data = dat_compt_10_12_adu, family = binomial, 
                                      control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between group and agency
compt_mod_10_12_adu_drop_Ac_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                         percent_comprehension + wob + gender_coded + group + 
                                         group:delta_EV + group:delta_uncertainty_level +
                                         (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                       data = dat_compt_10_12_adu, family = binomial, 
                                       control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop EV fixed
compt_mod_10_12_adu_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                           percent_comprehension + wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_10_12_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop uncertainty fixed
compt_mod_10_12_adu_drop_C_fixed <- glmer(info_choice ~ delta_EV + delta_agency +
                                          percent_comprehension + wob + gender_coded + group + 
                                          group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                          (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                        data = dat_compt_10_12_adu, family = binomial, 
                                        control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop agency fixed
compt_mod_10_12_adu_drop_Ac_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                           percent_comprehension + wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_10_12_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop comprehension fixed
compt_mod_10_12_adu_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                           wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_compt_10_12_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop wob fixed
compt_mod_10_12_adu_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                            percent_comprehension + gender_coded + group + 
                                            group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                            (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                          data = dat_compt_10_12_adu, family = binomial, 
                                          control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop gender fixed 
compt_mod_10_12_adu_drop_gender_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob + group + 
                                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                             data = dat_compt_10_12_adu, family = binomial, 
                                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop group fixed
compt_mod_10_12_adu_drop_group_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                           data = dat_compt_10_12_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)


# Model results
sjPlot::tab_model(compt_mod_child_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
sjPlot::tab_model(compt_mod_4_5_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
sjPlot::tab_model(compt_mod_6_7_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
sjPlot::tab_model(compt_mod_8_9_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE,show.ci=FALSE, show.se=TRUE)
sjPlot::tab_model(compt_mod_10_12_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)

## Test interactions
# Children versus adults
test_child_adu_EV_int <- anova(compt_mod_child_adu, compt_mod_child_adu_drop_Af_int, test="Chisq")
test_child_adu_uncertainty_int <- anova(compt_mod_child_adu, compt_mod_child_adu_drop_C_int, test="Chisq")
test_child_adu_agency_int <- anova(compt_mod_child_adu, compt_mod_child_adu_drop_Ac_int, test="Chisq")
test_child_adu_EV_fixed <- anova(compt_mod_child_adu, compt_mod_child_adu_drop_Af_fixed, test="Chisq")
test_child_adu_uncertainty_fixed <- anova(compt_mod_child_adu, compt_mod_child_adu_drop_C_fixed, test="Chisq")
test_child_adu_agency_fixed <- anova(compt_mod_child_adu, compt_mod_child_adu_drop_Ac_fixed, test="Chisq")
test_child_adu_pc_fixed <- anova(compt_mod_child_adu, compt_mod_child_adu_drop_pc_fixed, test="Chisq")
test_child_adu_wob_fixed <- anova(compt_mod_child_adu, compt_mod_child_adu_drop_wob_fixed, test="Chisq")
test_child_adu_gender_fixed <- anova(compt_mod_child_adu, compt_mod_child_adu_drop_gender_fixed, test="Chisq")
test_child_adu_group_fixed <- anova(compt_mod_child_adu, compt_mod_child_adu_drop_group_fixed, test="Chisq")

# 4-5 versus adults
test_4_5_adu_EV_int <- anova(compt_mod_4_5_adu, compt_mod_4_5_adu_drop_Af_int, test="Chisq")
test_4_5_adu_uncertainty_int <- anova(compt_mod_4_5_adu, compt_mod_4_5_adu_drop_C_int, test="Chisq")
test_4_5_adu_agency_int <- anova(compt_mod_4_5_adu, compt_mod_4_5_adu_drop_Ac_int, test="Chisq")
test_4_5_adu_EV_fixed <- anova(compt_mod_4_5_adu, compt_mod_4_5_adu_drop_Af_fixed, test="Chisq")
test_4_5_adu_uncertainty_fixed <- anova(compt_mod_4_5_adu, compt_mod_4_5_adu_drop_C_fixed, test="Chisq")
test_4_5_adu_agency_fixed <- anova(compt_mod_4_5_adu, compt_mod_4_5_adu_drop_Ac_fixed, test="Chisq")
test_4_5_adu_pc_fixed <- anova(compt_mod_4_5_adu, compt_mod_4_5_adu_drop_pc_fixed, test="Chisq")
test_4_5_adu_wob_fixed <- anova(compt_mod_4_5_adu, compt_mod_4_5_adu_drop_wob_fixed, test="Chisq")
test_4_5_adu_gender_fixed <- anova(compt_mod_4_5_adu, compt_mod_4_5_adu_drop_gender_fixed, test="Chisq")
test_4_5_adu_group_fixed <- anova(compt_mod_4_5_adu, compt_mod_4_5_adu_drop_group_fixed, test="Chisq")

# 6-7 versus adults
test_6_7_adu_EV_int <- anova(compt_mod_6_7_adu, compt_mod_6_7_adu_drop_Af_int, test="Chisq")
test_6_7_adu_uncertainty_int <- anova(compt_mod_6_7_adu, compt_mod_6_7_adu_drop_C_int, test="Chisq")
test_6_7_adu_agency_int <- anova(compt_mod_6_7_adu, compt_mod_6_7_adu_drop_Ac_int, test="Chisq")
test_6_7_adu_EV_fixed <- anova(compt_mod_6_7_adu, compt_mod_6_7_adu_drop_Af_fixed, test="Chisq")
test_6_7_adu_uncertainty_fixed <- anova(compt_mod_6_7_adu, compt_mod_6_7_adu_drop_C_fixed, test="Chisq")
test_6_7_adu_agency_fixed <- anova(compt_mod_6_7_adu, compt_mod_6_7_adu_drop_Ac_fixed, test="Chisq")
test_6_7_adu_pc_fixed <- anova(compt_mod_6_7_adu, compt_mod_6_7_adu_drop_pc_fixed, test="Chisq")
test_6_7_adu_wob_fixed <- anova(compt_mod_6_7_adu, compt_mod_6_7_adu_drop_wob_fixed, test="Chisq")
test_6_7_adu_gender_fixed <- anova(compt_mod_6_7_adu, compt_mod_6_7_adu_drop_gender_fixed, test="Chisq")
test_6_7_adu_group_fixed <- anova(compt_mod_6_7_adu, compt_mod_6_7_adu_drop_group_fixed, test="Chisq")

# 8-9 versus adults
test_8_9_adu_EV_int <- anova(compt_mod_8_9_adu, compt_mod_8_9_adu_drop_Af_int, test="Chisq")
test_8_9_adu_uncertainty_int <- anova(compt_mod_8_9_adu, compt_mod_8_9_adu_drop_C_int, test="Chisq")
test_8_9_adu_agency_int <- anova(compt_mod_8_9_adu, compt_mod_8_9_adu_drop_Ac_int, test="Chisq")
test_8_9_adu_EV_fixed <- anova(compt_mod_8_9_adu, compt_mod_8_9_adu_drop_Af_fixed, test="Chisq")
test_8_9_adu_uncertainty_fixed <- anova(compt_mod_8_9_adu, compt_mod_8_9_adu_drop_C_fixed, test="Chisq")
test_8_9_adu_agency_fixed <- anova(compt_mod_8_9_adu, compt_mod_8_9_adu_drop_Ac_fixed, test="Chisq")
test_8_9_adu_pc_fixed <- anova(compt_mod_8_9_adu, compt_mod_8_9_adu_drop_pc_fixed, test="Chisq")
test_8_9_adu_wob_fixed <- anova(compt_mod_8_9_adu, compt_mod_8_9_adu_drop_wob_fixed, test="Chisq")
test_8_9_adu_gender_fixed <- anova(compt_mod_8_9_adu, compt_mod_8_9_adu_drop_gender_fixed, test="Chisq")
test_8_9_adu_group_fixed <- anova(compt_mod_8_9_adu, compt_mod_8_9_adu_drop_group_fixed, test="Chisq")

# 10-12 versus adults
test_10_12_adu_EV_int <- anova(compt_mod_10_12_adu, compt_mod_10_12_adu_drop_Af_int, test="Chisq")
test_10_12_adu_uncertainty_int <- anova(compt_mod_10_12_adu, compt_mod_10_12_adu_drop_C_int, test="Chisq")
test_10_12_adu_agency_int <- anova(compt_mod_10_12_adu, compt_mod_10_12_adu_drop_Ac_int, test="Chisq")
test_10_12_adu_EV_fixed <- anova(compt_mod_10_12_adu, compt_mod_10_12_adu_drop_Af_fixed, test="Chisq")
test_10_12_adu_uncertainty_fixed <- anova(compt_mod_10_12_adu, compt_mod_10_12_adu_drop_C_fixed, test="Chisq")
test_10_12_adu_agency_fixed <- anova(compt_mod_10_12_adu, compt_mod_10_12_adu_drop_Ac_fixed, test="Chisq")
test_10_12_adu_pc_fixed <- anova(compt_mod_10_12_adu, compt_mod_10_12_adu_drop_pc_fixed, test="Chisq")
test_10_12_adu_wob_fixed <- anova(compt_mod_10_12_adu, compt_mod_10_12_adu_drop_wob_fixed, test="Chisq")
test_10_12_adu_gender_fixed <- anova(compt_mod_10_12_adu, compt_mod_10_12_adu_drop_gender_fixed, test="Chisq")
test_10_12_adu_group_fixed <- anova(compt_mod_10_12_adu, compt_mod_10_12_adu_drop_group_fixed, test="Chisq")


# Create csv file with tests of fixed effects
variable <- c("child_adu_EV",  "child_adu_uncertainty", "child_adu_agency", 
              "child_adu_comprehension", "child_adu_wob", "child_adu_gender", 
              "child_adu_group", "child_adu_group_by_EV", 
              "child_adu_group_by_uncertainty", "child_adu_group_by_agency",
              "4_5_adu_EV",  "4_5_adu_uncertainty", "4_5_adu_agency", 
              "4_5_adu_comprehension", "4_5_adu_wob", "4_5_adu_gender", 
              "4_5_adu_group", "4_5_adu_group_by_EV", 
              "4_5_adu_group_by_uncertainty", "4_5_adu_group_by_agency",
              "6_7_adu_EV",  "6_7_adu_uncertainty", "6_7_adu_agency", 
              "6_7_adu_comprehension", "6_7_adu_wob", "6_7_adu_gender", 
              "6_7_adu_group", "6_7_adu_group_by_EV", 
              "6_7_adu_group_by_uncertainty", "6_7_adu_group_by_agency",
              "8_9_adu_EV",  "8_9_adu_uncertainty", "8_9_adu_agency", 
              "8_9_adu_comprehension", "8_9_adu_wob", "8_9_adu_gender", 
              "8_9_adu_group", "8_9_adu_group_by_EV", 
              "8_9_adu_group_by_uncertainty", "8_9_adu_group_by_agency",
              "10_12_adu_EV",  "10_12_adu_uncertainty", "10_12_adu_agency", 
              "10_12_adu_comprehension", "10_12_adu_wob", "10_12_adu_gender", 
              "10_12_adu_group", "10_12_adu_group_by_EV", 
              "10_12_adu_group_by_uncertainty", "10_12_adu_group_by_agency")

# store coefficients
child_adu_compt_mod_full_coeffs <- summary(compt_mod_child_adu)$coefficients[, 1]
compt_mod_4_5_adu_full_coeffs <- summary(compt_mod_4_5_adu)$coefficients[, 1]
compt_mod_6_7_adu_full_coeffs <- summary(compt_mod_6_7_adu)$coefficients[, 1]
compt_mod_8_9_adu_full_coeffs <- summary(compt_mod_8_9_adu)$coefficients[, 1]
compt_mod_10_12_adu_full_coeffs <- summary(compt_mod_10_12_adu)$coefficients[, 1]

# store sems
child_adu_compt_mod_full_sems <- summary(compt_mod_child_adu)$coefficients[, 2]
compt_mod_4_5_adu_full_sems <- summary(compt_mod_4_5_adu)$coefficients[, 2]
compt_mod_6_7_adu_full_sems <- summary(compt_mod_6_7_adu)$coefficients[, 2]
compt_mod_8_9_adu_full_sems <- summary(compt_mod_8_9_adu)$coefficients[, 2]
compt_mod_10_12_adu_full_sems <- summary(compt_mod_10_12_adu)$coefficients[, 2]


estimate <- c(child_adu_compt_mod_full_coeffs["delta_EV"], 
              child_adu_compt_mod_full_coeffs["delta_uncertainty_level"],
              child_adu_compt_mod_full_coeffs["delta_agency"],
              child_adu_compt_mod_full_coeffs["percent_comprehension"],
              child_adu_compt_mod_full_coeffs["wob"],
              child_adu_compt_mod_full_coeffs["gender_coded1"],
              child_adu_compt_mod_full_coeffs["group"],
              child_adu_compt_mod_full_coeffs["delta_EV:group"],
              child_adu_compt_mod_full_coeffs["delta_uncertainty_level:group"],
              child_adu_compt_mod_full_coeffs["delta_agency:group"],
              compt_mod_4_5_adu_full_coeffs["delta_EV"], 
              compt_mod_4_5_adu_full_coeffs["delta_uncertainty_level"],
              compt_mod_4_5_adu_full_coeffs["delta_agency"],
              compt_mod_4_5_adu_full_coeffs["percent_comprehension"],
              compt_mod_4_5_adu_full_coeffs["wob"],
              compt_mod_4_5_adu_full_coeffs["gender_coded1"],
              compt_mod_4_5_adu_full_coeffs["group"],
              compt_mod_4_5_adu_full_coeffs["delta_EV:group"],
              compt_mod_4_5_adu_full_coeffs["delta_uncertainty_level:group"],
              compt_mod_4_5_adu_full_coeffs["delta_agency:group"],
              compt_mod_6_7_adu_full_coeffs["delta_EV"], 
              compt_mod_6_7_adu_full_coeffs["delta_uncertainty_level"],
              compt_mod_6_7_adu_full_coeffs["delta_agency"],
              compt_mod_6_7_adu_full_coeffs["percent_comprehension"],
              compt_mod_6_7_adu_full_coeffs["wob"],
              compt_mod_6_7_adu_full_coeffs["gender_coded1"],
              compt_mod_6_7_adu_full_coeffs["group"],
              compt_mod_6_7_adu_full_coeffs["delta_EV:group"],
              compt_mod_6_7_adu_full_coeffs["delta_uncertainty_level:group"],
              compt_mod_6_7_adu_full_coeffs["delta_agency:group"],
              compt_mod_8_9_adu_full_coeffs["delta_EV"], 
              compt_mod_8_9_adu_full_coeffs["delta_uncertainty_level"],
              compt_mod_8_9_adu_full_coeffs["delta_agency"],
              compt_mod_8_9_adu_full_coeffs["percent_comprehension"],
              compt_mod_8_9_adu_full_coeffs["wob"],
              compt_mod_8_9_adu_full_coeffs["gender_coded1"],
              compt_mod_8_9_adu_full_coeffs["group"],
              compt_mod_8_9_adu_full_coeffs["delta_EV:group"],
              compt_mod_8_9_adu_full_coeffs["delta_uncertainty_level:group"],
              compt_mod_8_9_adu_full_coeffs["delta_agency:group"],
              compt_mod_10_12_adu_full_coeffs["delta_EV"], 
              compt_mod_10_12_adu_full_coeffs["delta_uncertainty_level"],
              compt_mod_10_12_adu_full_coeffs["delta_agency"],
              compt_mod_10_12_adu_full_coeffs["percent_comprehension"],
              compt_mod_10_12_adu_full_coeffs["wob"],
              compt_mod_10_12_adu_full_coeffs["gender_coded1"],
              compt_mod_10_12_adu_full_coeffs["group"],
              compt_mod_10_12_adu_full_coeffs["delta_EV:group"],
              compt_mod_10_12_adu_full_coeffs["delta_uncertainty_level:group"],
              compt_mod_10_12_adu_full_coeffs["delta_agency:group"])

sem <- c(child_adu_compt_mod_full_sems["delta_EV"], 
              child_adu_compt_mod_full_sems["delta_uncertainty_level"],
              child_adu_compt_mod_full_sems["delta_agency"],
              child_adu_compt_mod_full_sems["percent_comprehension"],
              child_adu_compt_mod_full_sems["wob"],
              child_adu_compt_mod_full_sems["gender_coded1"],
              child_adu_compt_mod_full_sems["group"],
              child_adu_compt_mod_full_sems["delta_EV:group"],
              child_adu_compt_mod_full_sems["delta_uncertainty_level:group"],
              child_adu_compt_mod_full_sems["delta_agency:group"],
              compt_mod_4_5_adu_full_sems["delta_EV"], 
              compt_mod_4_5_adu_full_sems["delta_uncertainty_level"],
              compt_mod_4_5_adu_full_sems["delta_agency"],
              compt_mod_4_5_adu_full_sems["percent_comprehension"],
              compt_mod_4_5_adu_full_sems["wob"],
              compt_mod_4_5_adu_full_sems["gender_coded1"],
              compt_mod_4_5_adu_full_sems["group"],
              compt_mod_4_5_adu_full_sems["delta_EV:group"],
              compt_mod_4_5_adu_full_sems["delta_uncertainty_level:group"],
              compt_mod_4_5_adu_full_sems["delta_agency:group"],
              compt_mod_6_7_adu_full_sems["delta_EV"], 
              compt_mod_6_7_adu_full_sems["delta_uncertainty_level"],
              compt_mod_6_7_adu_full_sems["delta_agency"],
              compt_mod_6_7_adu_full_sems["percent_comprehension"],
              compt_mod_6_7_adu_full_sems["wob"],
              compt_mod_6_7_adu_full_sems["gender_coded1"],
              compt_mod_6_7_adu_full_sems["group"],
              compt_mod_6_7_adu_full_sems["delta_EV:group"],
              compt_mod_6_7_adu_full_sems["delta_uncertainty_level:group"],
              compt_mod_6_7_adu_full_sems["delta_agency:group"],
              compt_mod_8_9_adu_full_sems["delta_EV"], 
              compt_mod_8_9_adu_full_sems["delta_uncertainty_level"],
              compt_mod_8_9_adu_full_sems["delta_agency"],
              compt_mod_8_9_adu_full_sems["percent_comprehension"],
              compt_mod_8_9_adu_full_sems["wob"],
              compt_mod_8_9_adu_full_sems["gender_coded1"],
              compt_mod_8_9_adu_full_sems["group"],
              compt_mod_8_9_adu_full_sems["delta_EV:group"],
              compt_mod_8_9_adu_full_sems["delta_uncertainty_level:group"],
              compt_mod_8_9_adu_full_sems["delta_agency:group"],
              compt_mod_10_12_adu_full_sems["delta_EV"], 
              compt_mod_10_12_adu_full_sems["delta_uncertainty_level"],
              compt_mod_10_12_adu_full_sems["delta_agency"],
              compt_mod_10_12_adu_full_sems["percent_comprehension"],
              compt_mod_10_12_adu_full_sems["wob"],
              compt_mod_10_12_adu_full_sems["gender_coded1"],
              compt_mod_10_12_adu_full_sems["group"],
              compt_mod_10_12_adu_full_sems["delta_EV:group"],
              compt_mod_10_12_adu_full_sems["delta_uncertainty_level:group"],
              compt_mod_10_12_adu_full_sems["delta_agency:group"])

chisq <- c(test_child_adu_EV_fixed$Chisq[2], test_child_adu_uncertainty_fixed$Chisq[2],
           test_child_adu_agency_fixed$Chisq[2], test_child_adu_pc_fixed$Chisq[2],
           test_child_adu_wob_fixed$Chisq[2], test_child_adu_gender_fixed$Chisq[2], 
           test_child_adu_group_fixed$Chisq[2], test_child_adu_EV_int$Chisq[2], 
           test_child_adu_uncertainty_int$Chisq[2], test_child_adu_agency_int$Chisq[2],
           test_4_5_adu_EV_fixed$Chisq[2], test_4_5_adu_uncertainty_fixed$Chisq[2],
           test_4_5_adu_agency_fixed$Chisq[2], test_4_5_adu_pc_fixed$Chisq[2],
           test_4_5_adu_wob_fixed$Chisq[2], test_4_5_adu_gender_fixed$Chisq[2], 
           test_4_5_adu_group_fixed$Chisq[2], test_4_5_adu_EV_int$Chisq[2], 
           test_4_5_adu_uncertainty_int$Chisq[2], test_4_5_adu_agency_int$Chisq[2],
           test_6_7_adu_EV_fixed$Chisq[2], test_6_7_adu_uncertainty_fixed$Chisq[2],
           test_6_7_adu_agency_fixed$Chisq[2], test_6_7_adu_pc_fixed$Chisq[2],
           test_6_7_adu_wob_fixed$Chisq[2], test_6_7_adu_gender_fixed$Chisq[2], 
           test_6_7_adu_group_fixed$Chisq[2], test_6_7_adu_EV_int$Chisq[2], 
           test_6_7_adu_uncertainty_int$Chisq[2], test_6_7_adu_agency_int$Chisq[2],
           test_8_9_adu_EV_fixed$Chisq[2], test_8_9_adu_uncertainty_fixed$Chisq[2],
           test_8_9_adu_agency_fixed$Chisq[2], test_8_9_adu_pc_fixed$Chisq[2],
           test_8_9_adu_wob_fixed$Chisq[2], test_8_9_adu_gender_fixed$Chisq[2], 
           test_8_9_adu_group_fixed$Chisq[2], test_8_9_adu_EV_int$Chisq[2], 
           test_8_9_adu_uncertainty_int$Chisq[2], test_8_9_adu_agency_int$Chisq[2],
           test_10_12_adu_EV_fixed$Chisq[2], test_10_12_adu_uncertainty_fixed$Chisq[2],
           test_10_12_adu_agency_fixed$Chisq[2], test_10_12_adu_pc_fixed$Chisq[2],
           test_10_12_adu_wob_fixed$Chisq[2], test_10_12_adu_gender_fixed$Chisq[2], 
           test_10_12_adu_group_fixed$Chisq[2], test_10_12_adu_EV_int$Chisq[2], 
           test_10_12_adu_uncertainty_int$Chisq[2], test_10_12_adu_agency_int$Chisq[2])

df <- c(test_child_adu_EV_fixed$Df[2], test_child_adu_uncertainty_fixed$Df[2],
           test_child_adu_agency_fixed$Df[2], test_child_adu_pc_fixed$Df[2],
           test_child_adu_wob_fixed$Df[2], test_child_adu_gender_fixed$Df[2], 
           test_child_adu_group_fixed$Df[2], test_child_adu_EV_int$Df[2], 
           test_child_adu_uncertainty_int$Df[2], test_child_adu_agency_int$Df[2],
           test_4_5_adu_EV_fixed$Df[2], test_4_5_adu_uncertainty_fixed$Df[2],
           test_4_5_adu_agency_fixed$Df[2], test_4_5_adu_pc_fixed$Df[2],
           test_4_5_adu_wob_fixed$Df[2], test_4_5_adu_gender_fixed$Df[2], 
           test_4_5_adu_group_fixed$Df[2], test_4_5_adu_EV_int$Df[2], 
           test_4_5_adu_uncertainty_int$Df[2], test_4_5_adu_agency_int$Df[2],
           test_6_7_adu_EV_fixed$Df[2], test_6_7_adu_uncertainty_fixed$Df[2],
           test_6_7_adu_agency_fixed$Df[2], test_6_7_adu_pc_fixed$Df[2],
           test_6_7_adu_wob_fixed$Df[2], test_6_7_adu_gender_fixed$Df[2], 
           test_6_7_adu_group_fixed$Df[2], test_6_7_adu_EV_int$Df[2], 
           test_6_7_adu_uncertainty_int$Df[2], test_6_7_adu_agency_int$Df[2],
           test_8_9_adu_EV_fixed$Df[2], test_8_9_adu_uncertainty_fixed$Df[2],
           test_8_9_adu_agency_fixed$Df[2], test_8_9_adu_pc_fixed$Df[2],
           test_8_9_adu_wob_fixed$Df[2], test_8_9_adu_gender_fixed$Df[2], 
           test_8_9_adu_group_fixed$Df[2], test_8_9_adu_EV_int$Df[2], 
           test_8_9_adu_uncertainty_int$Df[2], test_8_9_adu_agency_int$Df[2],
           test_10_12_adu_EV_fixed$Df[2], test_10_12_adu_uncertainty_fixed$Df[2],
           test_10_12_adu_agency_fixed$Df[2], test_10_12_adu_pc_fixed$Df[2],
           test_10_12_adu_wob_fixed$Df[2], test_10_12_adu_gender_fixed$Df[2], 
           test_10_12_adu_group_fixed$Df[2], test_10_12_adu_EV_int$Df[2], 
           test_10_12_adu_uncertainty_int$Df[2], test_10_12_adu_agency_int$Df[2])

p_value <- c(test_child_adu_EV_fixed$Pr[2], test_child_adu_uncertainty_fixed$Pr[2],
           test_child_adu_agency_fixed$Pr[2], test_child_adu_pc_fixed$Pr[2],
           test_child_adu_wob_fixed$Pr[2], test_child_adu_gender_fixed$Pr[2], 
           test_child_adu_group_fixed$Pr[2], test_child_adu_EV_int$Pr[2], 
           test_child_adu_uncertainty_int$Pr[2], test_child_adu_agency_int$Pr[2],
           test_4_5_adu_EV_fixed$Pr[2], test_4_5_adu_uncertainty_fixed$Pr[2],
           test_4_5_adu_agency_fixed$Pr[2], test_4_5_adu_pc_fixed$Pr[2],
           test_4_5_adu_wob_fixed$Pr[2], test_4_5_adu_gender_fixed$Pr[2], 
           test_4_5_adu_group_fixed$Pr[2], test_4_5_adu_EV_int$Pr[2], 
           test_4_5_adu_uncertainty_int$Pr[2], test_4_5_adu_agency_int$Pr[2],
           test_6_7_adu_EV_fixed$Pr[2], test_6_7_adu_uncertainty_fixed$Pr[2],
           test_6_7_adu_agency_fixed$Pr[2], test_6_7_adu_pc_fixed$Pr[2],
           test_6_7_adu_wob_fixed$Pr[2], test_6_7_adu_gender_fixed$Pr[2], 
           test_6_7_adu_group_fixed$Pr[2], test_6_7_adu_EV_int$Pr[2], 
           test_6_7_adu_uncertainty_int$Pr[2], test_6_7_adu_agency_int$Pr[2],
           test_8_9_adu_EV_fixed$Pr[2], test_8_9_adu_uncertainty_fixed$Pr[2],
           test_8_9_adu_agency_fixed$Pr[2], test_8_9_adu_pc_fixed$Pr[2],
           test_8_9_adu_wob_fixed$Pr[2], test_8_9_adu_gender_fixed$Pr[2], 
           test_8_9_adu_group_fixed$Pr[2], test_8_9_adu_EV_int$Pr[2], 
           test_8_9_adu_uncertainty_int$Pr[2], test_8_9_adu_agency_int$Pr[2],
           test_10_12_adu_EV_fixed$Pr[2], test_10_12_adu_uncertainty_fixed$Pr[2],
           test_10_12_adu_agency_fixed$Pr[2], test_10_12_adu_pc_fixed$Pr[2],
           test_10_12_adu_wob_fixed$Pr[2], test_10_12_adu_gender_fixed$Pr[2], 
           test_10_12_adu_group_fixed$Pr[2], test_10_12_adu_EV_int$Pr[2], 
           test_10_12_adu_uncertainty_int$Pr[2], test_10_12_adu_agency_int$Pr[2])


compt_mods_fixed_effects_child_vs_adu <- data.frame(variable, estimate, sem, chisq, df, p_value)

# Save fixed effects as a file
write.csv(compt_mods_fixed_effects_child_vs_adu, sprintf("%s/compt_mods_fixed_effects_child_vs_adu.csv", wd), row.names = FALSE)



#### Fishing choices ####
# t-tests
# based on non-decoy items
t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_w_info_non_Z_1, 
            subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_no_info_non_Z_1,
       var.equal=TRUE)

t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_w_info_non_Z_1, 
       subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_no_info_non_Z_1,
       var.equal=TRUE)

t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_w_info_non_Z_1, 
       subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_no_info_non_Z_1,
       var.equal=TRUE)

t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_w_info_non_Z_1, 
       subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_no_info_non_Z_1,
       var.equal=TRUE)

t.test(subset(adu_dat_compt, trial_number == 1)$fishing_w_info_non_Z_1, 
       subset(adu_dat_compt, trial_number == 1)$fishing_no_info_non_Z_1,
       var.equal=TRUE)

# based on what participants see
t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_w_info_non_Z_2, 
       subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_no_info_non_Z_2,
       var.equal=TRUE)

t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_w_info_non_Z_2, 
       subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_no_info_non_Z_2,
       var.equal=TRUE)

t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_w_info_non_Z_2, 
       subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_no_info_non_Z_2,
       var.equal=TRUE)

t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_w_info_non_Z_2, 
       subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_no_info_non_Z_2,
       var.equal=TRUE)

t.test(subset(adu_dat_compt, trial_number == 1)$fishing_w_info_non_Z_2, 
       subset(adu_dat_compt, trial_number == 1)$fishing_no_info_non_Z_2,
       var.equal=TRUE)

# Overall fishing choices compared to chance and children vs adults
# based on what participants see
t.test(subset(dat_compt, trial_number == 1)$fishing_non_Z_2, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(dat_compt, trial_number == 1)$fishing_non_Z_2)
t.test(subset(adu_dat_compt, trial_number == 1)$fishing_non_Z_2, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(adu_dat_compt, trial_number == 1)$fishing_non_Z_2)
t.test(subset(dat_compt, trial_number == 1)$fishing_non_Z_2, 
       subset(adu_dat_compt, trial_number == 1)$fishing_non_Z_2,var.equal=TRUE)

t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_non_Z_2, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_non_Z_2)
t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_non_Z_2, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_non_Z_2)
t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_non_Z_2, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_non_Z_2)
t.test(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_non_Z_2, 
       mu=0.5, var.equal=TRUE)
sdamr::sample_sd(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_non_Z_2)

# Plot average scores by age group
# Function to get the standard error from the mean
sem <- function(x) sdamr::sample_sd(x)/sqrt(length(x))
# Create dataset for fishing choices
age_group <- factor(c("4-5", "6-7", "8-9", "10-12", "Adults"), levels=c("4-5", "6-7", "8-9", "10-12", "Adults"))
x_points <- c(0.1, 0.2, 0.3, 0.4, 0.5)
x_labels <- c("4-5", "6-7", "8-9", "10-12", "Adults")

# based on non-decoy items
fishing_w_info_1 <- c(
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_w_info_non_Z_1),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_w_info_non_Z_1),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_w_info_non_Z_1),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_w_info_non_Z_1),
  mean(subset(adu_dat_compt, trial_number == 1)$fishing_w_info_non_Z_1)
)
fishing_w_info_sem_1 <- c(
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_w_info_non_Z_1),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_w_info_non_Z_1),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_w_info_non_Z_1),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_w_info_non_Z_1),
  sem(subset(adu_dat_compt, trial_number == 1)$fishing_w_info_non_Z_1)
)
fishing_no_info_1 <- c(
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_no_info_non_Z_1),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_no_info_non_Z_1),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_no_info_non_Z_1),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_no_info_non_Z_1),
  mean(subset(adu_dat_compt, trial_number == 1)$fishing_no_info_non_Z_1)
)
fishing_no_info_sem_1 <- c(
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_no_info_non_Z_1),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_no_info_non_Z_1),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_no_info_non_Z_1),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_no_info_non_Z_1),
  sem(subset(adu_dat_compt, trial_number == 1)$fishing_no_info_non_Z_1)
)
fishing_1 <- c(
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_non_Z_1),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_non_Z_1),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_non_Z_1),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_non_Z_1),
  mean(subset(adu_dat_compt, trial_number == 1)$fishing_non_Z_1)
)
fishing_sem_1 <- c(
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_non_Z_1),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_non_Z_1),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_non_Z_1),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_non_Z_1),
  sem(subset(adu_dat_compt, trial_number == 1)$fishing_non_Z_1)
)

# based on what participants see
fishing_w_info_2 <- c(
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_w_info_non_Z_2),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_w_info_non_Z_2),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_w_info_non_Z_2),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_w_info_non_Z_2),
  mean(subset(adu_dat_compt, trial_number == 1)$fishing_w_info_non_Z_2)
)
fishing_w_info_sem_2 <- c(
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_w_info_non_Z_2),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_w_info_non_Z_2),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_w_info_non_Z_2),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_w_info_non_Z_2),
  sem(subset(adu_dat_compt, trial_number == 1)$fishing_w_info_non_Z_2)
)
fishing_no_info_2 <- c(
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_no_info_non_Z_2),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_no_info_non_Z_2),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_no_info_non_Z_2),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_no_info_non_Z_2),
  mean(subset(adu_dat_compt, trial_number == 1)$fishing_no_info_non_Z_2)
)
fishing_no_info_sem_2 <- c(
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_no_info_non_Z_2),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_no_info_non_Z_2),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_no_info_non_Z_2),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_no_info_non_Z_2),
  sem(subset(adu_dat_compt, trial_number == 1)$fishing_no_info_non_Z_2)
)
fishing_2 <- c(
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_non_Z_2),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_non_Z_2),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_non_Z_2),
  mean(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_non_Z_2),
  mean(subset(adu_dat_compt, trial_number == 1)$fishing_non_Z_2)
)
fishing_sem_2 <- c(
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 0)$fishing_non_Z_2),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 1)$fishing_non_Z_2),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 2)$fishing_non_Z_2),
  sem(subset(dat_compt, trial_number == 1 & age_group_coded == 3)$fishing_non_Z_2),
  sem(subset(adu_dat_compt, trial_number == 1)$fishing_non_Z_2)
)

fishing_choices_dat_by_group <- data.frame(x_points, x_labels, fishing_w_info_1,
                                           fishing_w_info_sem_1, fishing_no_info_1,
                                           fishing_1, fishing_sem_2,
                                           fishing_no_info_sem_2, fishing_w_info_2,
                                           fishing_w_info_sem_2, fishing_no_info_2,
                                           fishing_no_info_sem_2, 
                                           fishing_2, fishing_sem_2)

# based on non-decoy items
fishing_choices_1_by_group_plot <- ggplot(fishing_choices_dat_by_group) + 
  geom_line(aes(y=fishing_w_info_1, x=x_points), size=1, color="orange1")+
  geom_point(aes(y=fishing_w_info_1, x=x_points), color="orange1", size = 4)+
  geom_ribbon(aes(ymin=fishing_w_info_1-fishing_w_info_sem_1, 
                  ymax=fishing_w_info_1+fishing_w_info_sem_1, 
                  x=x_points), fill="orange1" , alpha = 0.3)+
  geom_line(aes(y=fishing_no_info_1, x=x_points), color="purple3", size=1)+
  geom_point(aes(y=fishing_no_info_1, x=x_points), color="purple3", size = 4)+
  geom_ribbon(aes(ymin=fishing_no_info_1-fishing_no_info_sem_1, 
                  ymax=fishing_no_info_1+fishing_no_info_sem_1, 
                  x=x_points), fill="purple3", alpha = 0.3)+
  theme_classic()+
  theme(text = element_text(size=16), legend.position="none") +
  scale_x_continuous(labels=fishing_choices_dat_by_group$x_labels) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  geom_hline(yintercept=0.5, color="#AAAAAA") +
  labs(x="Age group", y="Correct fishing choices\n(based on non decoy items)")

fishing_choices_1_by_group_plot

fishing_choices_1_by_group_overall_plot <- ggplot(fishing_choices_dat_by_group) + 
  geom_line(aes(y=fishing_1, x=x_points), size=1, color="darkgrey")+
  geom_point(aes(y=fishing_1, x=x_points), color="darkgrey", size = 4)+
  geom_ribbon(aes(ymin=fishing_1-fishing_sem_1, 
                  ymax=fishing_1+fishing_sem_1, 
                  x=x_points), fill="darkgrey" , alpha = 0.3)+
  theme_classic()+
  theme(text = element_text(size=16), legend.position="none") +
  scale_x_continuous(labels=fishing_choices_dat_by_group$x_labels) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  geom_hline(yintercept=0.5, color="#AAAAAA") +
  labs(x="Age group", y="Correct fishing choices\n(based on non decoy items)")

fishing_choices_1_by_group_overall_plot

# based on what participants see
fishing_choices_2_by_group_plot <- ggplot(fishing_choices_dat_by_group) + 
  geom_line(aes(y=fishing_w_info_2, x=x_points), size=1, color="orange1")+
  geom_point(aes(y=fishing_w_info_2, x=x_points), color="orange1", size = 4)+
  geom_ribbon(aes(ymin=fishing_w_info_2-fishing_w_info_sem_2, 
                  ymax=fishing_w_info_2+fishing_w_info_sem_2, 
                  x=x_points), fill="orange1" , alpha = 0.3)+
  geom_line(aes(y=fishing_no_info_2, x=x_points), color="purple3", size=1)+
  geom_point(aes(y=fishing_no_info_2, x=x_points), color="purple3", size = 4)+
  geom_ribbon(aes(ymin=fishing_no_info_2-fishing_no_info_sem_2, 
                  ymax=fishing_no_info_2+fishing_no_info_sem_2, 
                  x=x_points), fill="purple3", alpha = 0.3)+
  theme_classic()+
  theme(text = element_text(size=16), legend.position="none") +
  scale_x_continuous(labels=fishing_choices_dat_by_group$x_labels) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  geom_hline(yintercept=0.5, color="#AAAAAA") +
  labs(x="Age group", y="Correct fishing choices\n(based on what participants see)")

fishing_choices_2_by_group_plot

fishing_choices_2_by_group_overall_plot <- ggplot(fishing_choices_dat_by_group) + 
  geom_line(aes(y=fishing_2, x=x_points), size=1, color="darkgrey")+
  geom_point(aes(y=fishing_2, x=x_points), color="darkgrey", size = 4)+
  geom_ribbon(aes(ymin=fishing_2-fishing_sem_2, 
                  ymax=fishing_2+fishing_sem_2, 
                  x=x_points), fill="darkgrey" , alpha = 0.3)+
  theme_classic()+
  theme(text = element_text(size=16), legend.position="none") +
  scale_x_continuous(labels=fishing_choices_dat_by_group$x_labels) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  geom_hline(yintercept=0.5) +
  labs(x="Age group", y="Correct fishing choices\n(based on what participants see)")

fishing_choices_2_by_group_overall_plot

