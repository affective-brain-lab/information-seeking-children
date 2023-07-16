# AUTHOR: GAIA MOLINARO
# Data analysis for the main experiments
# in Molinaro, Cogliati Dezza, Buehler, Moutsiana & Sharot (in prep.)

# CLEAR WORKSPACE
rm(list = ls())

# SET WORKING DIRECTORY
wd <- "insert_your_wd_here"
setwd(wd)

packages <- c("ggplot2", "GGally", "ggpubr", "psych", "nls2", "RColorBrewer",
              "grid", "optimx", "ggcorrplot", "tidyverse", "glue", "gridExtra", 
              "lme4", "table1", "BayesFactor", "modelbased", "quest",
              "glmnet", "glmnetUtils", "broom", "MuMIn", "ggtext")
pacman::p_load(packages, character.only = TRUE)
#devtools::install_github("gaiamolinaro/interactions")
library(interactions)

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

# get omega (chi-square test effect size)
get_omega <- function(chi_sq, N) {
  return (sqrt(chi_sq/N))
}



#### EXPERIMENTS 1 and 2 ####
#### Load and prep the dataset #### 
rescale_variables <- TRUE  # whether to rescale the numeric variables
center_variables <- TRUE # whether to center numeric variables (except deltas)
rescale_age <- TRUE  # whether to rescale age
center_age <- TRUE # whether to center age

data_batch <- 1  # 1 or 2
dat <- read.csv(sprintf("data_experiment_%s.csv", data_batch))
dat <- subset(dat, !(condition %in% c("catch_1", "catch_2")))
dat <- subset(dat, catch_trials_score == 100)

# Create age groups
dat$age_group <- "None"
dat$age_group[dat$age_in_years  %in% c(4, 5)] <- "4-5"
dat$age_group[dat$age_in_years  %in% c(6, 7)] <- "6-7"
dat$age_group[dat$age_in_years  %in% c(8, 9)] <- "8-9"
dat$age_group[dat$age_in_years  %in% c(10, 11, 12)] <- "10-12"
dat$age_group <- factor(dat$age_group, levels=c("4-5", "6-7", "8-9","10-12"))

# Turn variables into appropriate data types
# Use scale() for standardized values
dat$delta_EV_non_Z <- as.numeric(dat$delta_EV)
dat$delta_uncertainty_level_non_Z <- as.numeric(dat$delta_uncertainty_level)
dat$delta_agency_non_Z <- as.numeric(dat$delta_agency)
dat$delta_SD <- as.numeric(dat$delta_SD)
dat$age_in_years_non_Z <- as.numeric(dat$age_in_years)
dat$age_in_years <- scale(as.numeric(dat$age_in_years), center=center_age, scale=rescale_age)
dat$age_in_months_non_Z <- as.numeric(dat$age_in_months)
dat$age_in_months_log <- scale(log(as.numeric(dat$age_in_months_non_Z)), center=center_age, scale=rescale_age)
dat$age_in_months <- scale(as.numeric(dat$age_in_months), center=center_age, scale=rescale_age)
dat$age_group_coded <- as.factor(dat$age_group_coded)
dat$RT_info_choice_non_Z <- as.numeric(dat$RT_info_choice)
dat$RT_info_choice_log <- scale(log(as.numeric(dat$RT_info_choice_non_Z)), center=center_variables, scale=rescale_variables)
dat$RT_info_choice <- scale(as.numeric(dat$RT_info_choice), center=center_variables, scale=rescale_variables)
dat$subject_ID <- as.factor(dat$subject_ID)
dat$gender_coded[dat$gender_coded==3] <- 1 # 1 = male, 2 = female, 3 = other
dat$gender_coded <- as.factor(dat$gender_coded)
contrasts(dat$gender_coded) <- contr.helmert(2)

# previous information-seeking choice
dat$prev_info_choice <- as.factor(quest::shifts_by(data=dat, vrb.nm="info_choice", grp.nm = "subject_ID", n=-1L)[[1]])
dat$info_choice <- as.factor(dat$info_choice) # 0 = left, 1 = right
dat$percent_comprehension_non_Z <- as.numeric(dat$percent_comprehension)
dat$percent_comprehension <- scale(as.numeric(dat$percent_comprehension), center=center_variables, scale=rescale_variables)
dat$wob_non_Z <- as.numeric(dat$wob)

dat$wob <- scale(as.numeric(dat$wob), center=center_variables, scale=rescale_variables)
dat$fishing_rewardL <- as.numeric(dat$fishing_rewardL)
dat$fishing_rewardR <- as.numeric(dat$fishing_rewardR)
dat$reward_external <- as.numeric(dat$fishing_rewardL + dat$fishing_rewardR)
dat$fishing_non_Z <- as.numeric(dat$prop_correct_fishing)
dat$fishing <- scale(as.numeric(dat$prop_correct_fishing), center=center_variables, scale=rescale_variables)

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
  dat <- rescale_deltas(dat)  
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

check_delta_signs(dat)

# Create subsets of data for each age group
dat_4_5 <- subset(dat, age_group == "4-5")
dat_6_7 <- subset(dat, age_group == "6-7")
dat_8_9 <- subset(dat, age_group == "8-9")
dat_10_12 <- subset(dat, age_group == "10-12")
# Re-Z-score and rescale variables
if (rescale_variables) {
  dat_4_5$percent_comprehension <- scale(dat_4_5$percent_comprehension_non_Z, center=center_variables, scale=rescale_variables)
  dat_6_7$percent_comprehension <- scale(dat_6_7$percent_comprehension_non_Z, center=center_variables, scale=rescale_variables)
  dat_8_9$percent_comprehension <- scale(dat_8_9$percent_comprehension_non_Z, center=center_variables, scale=rescale_variables)
  dat_10_12$percent_comprehension <- scale(dat_10_12$percent_comprehension_non_Z, center=center_variables, scale=rescale_variables)
  dat_4_5$wob <- scale(dat_4_5$wob_non_Z, center=center_variables, scale=rescale_variables)
  dat_6_7$wob <- scale(dat_6_7$wob_non_Z, center=center_variables, scale=rescale_variables)
  dat_8_9$wob <- scale(dat_8_9$wob_non_Z, center=center_variables, scale=rescale_variables)
  dat_10_12$wob <- scale(dat_10_12$wob_non_Z, center=center_variables, scale=rescale_variables)
  dat_4_5 <- rescale_deltas(dat_4_5)
  dat_6_7 <- rescale_deltas(dat_6_7)
  dat_10_12 <- rescale_deltas(dat_10_12) 
}

#### Demographics ####
summary(subset(dat, trial_number==1)$age_in_years_non_Z)
sdamr::sample_sd(subset(dat, trial_number==1)$age_in_years_non_Z)
psych::describeBy(subset(dat, trial_number==1)$age_in_years_non_Z, 
                  subset(dat, trial_number==1)$age_group)
table(subset(dat, trial_number==1)$age_group)
table(subset(dat, trial_number==1)$age_in_years_non_Z)
table(subset(dat, trial_number==1)$gender)
table(subset(dat, trial_number==1)$gender, 
      subset(dat, trial_number==1)$age_group)
table1::table1(~ gender + age_in_years_non_Z + percent_comprehension_non_Z + 
                 wob_non_Z + fishing_non_Z | age_group, 
               data=subset(dat, dat$trial_number == 1))

#### Comprehension, wob, and fishing scores ####
##### By age group #####
## t-tests
## t-tests
t.test(subset(dat, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(adu_dat, trial_number == 1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_4_5, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_6_7, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_8_9, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(dat_10_12, trial_number ==1)$wob_non_Z, mu=0.5, var.equal=TRUE)

# normality tests
shapiro.test(subset(dat_4_5, trial_number ==1)$wob_non_Z)
shapiro.test(subset(dat_6_7, trial_number ==1)$wob_non_Z)
shapiro.test(subset(dat_8_9, trial_number ==1)$wob_non_Z)
shapiro.test(subset(dat_10_12, trial_number ==1)$wob_non_Z)

# non-parametric tests
wilcox.test(subset(dat_4_5, trial_number ==1)$wob_non_Z, mu=0.5, exact=FALSE)
wilcox.test(subset(dat_6_7, trial_number ==1)$wob_non_Z, mu=0.5, exact=FALSE)
wilcox.test(subset(dat_8_9, trial_number ==1)$wob_non_Z, mu=0.5, exact=FALSE)
wilcox.test(subset(dat_10_12, trial_number ==1)$wob_non_Z, mu=0.5, exact=FALSE)

subset(dat_4_5, trial_number ==1) %>% rstatix::wilcox_effsize(wob_non_Z ~ 1, mu = 0)
subset(dat_6_7, trial_number ==1) %>% rstatix::wilcox_effsize(wob_non_Z ~ 1, mu = 0)
subset(dat_8_9, trial_number ==1) %>% rstatix::wilcox_effsize(wob_non_Z ~ 1, mu = 0)
subset(dat_10_12, trial_number ==1) %>% rstatix::wilcox_effsize(wob_non_Z ~ 1, mu = 0)

# Descriptive stats by age group
psych::describeBy(subset(dat, trial_number==1)$percent_comprehension_non_Z, 
                  subset(dat, trial_number==1)$age_group)
psych::describeBy(subset(dat, trial_number==1)$wob_non_Z, 
                  subset(dat, trial_number==1)$age_group)
psych::describeBy(subset(dat, trial_number==1)$fishing_non_Z, 
                  subset(dat, trial_number==1)$age_group)

## Create long form data sets
# Create separate data sets with variables for children and adults' pc and wob
dat_pc_wob <- dplyr::select(subset(dat, trial_number==1), 
                                  all_of(c("subject_ID", "age_group_coded", 
                                           "percent_comprehension_non_Z", 
                                           "wob_non_Z")))
dat_pc_wob$percent_comprehension_child <- dat_pc_wob$"percent_comprehension_non_Z"
dat_pc_wob$wob_child <- dat_pc_wob$"wob_non_Z"

adu_dat_pc_wob <- dplyr::select(subset(adu_dat, trial_number==1), 
                                      all_of(c("subject_ID", "age_group_coded", 
                                               "percent_comprehension_non_Z", 
                                               "wob_non_Z")))

# Merge datasets

# Turn into long form 
dat_pc_wob <- tidyr::gather(dat_pc_wob, key = "variable", value = "score", percent_comprehension_non_Z:wob_non_Z, factor_key=TRUE)
dat_pc_wob$color <- "#AAAAAA"

# Create a new variable to distinguish between each group of children
dat_pc_wob$variable2 <- "None"
dat_pc_wob$variable2[dat_pc_wob$age_group_coded == 0
                                     & dat_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_4_5"
dat_pc_wob$variable2[dat_pc_wob$age_group_coded == 0
                                     & dat_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_4_5"
dat_pc_wob$variable2[dat_pc_wob$age_group_coded == 1
                                     & dat_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_6_7"
dat_pc_wob$variable2[dat_pc_wob$age_group_coded == 1
                                     & dat_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_6_7"
dat_pc_wob$variable2[dat_pc_wob$age_group_coded == 2
                                     & dat_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_8_9"
dat_pc_wob$variable2[dat_pc_wob$age_group_coded == 2
                                     & dat_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_8_9"
dat_pc_wob$variable2[dat_pc_wob$age_group_coded == 3
                                     & dat_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_10_12"
dat_pc_wob$variable2[dat_pc_wob$age_group_coded == 3
                                     & dat_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_10_12"
dat_pc_wob$variable2 <- factor(dat_pc_wob$variable2, 
                                               levels=c("percent_comprehension_4_5",
                                                        "percent_comprehension_6_7",
                                                        "percent_comprehension_8_9",
                                                        "percent_comprehension_10_12",
                                                        "wob_4_5",
                                                        "wob_6_7",
                                                        "wob_8_9",
                                                        "wob_10_12"))

table(dat_pc_wob$variable2)

## Plots
x_labels <- c("4-5", "6-7", "8-9", "10-12")

# By group, without adults
pc_by_group <- ggplot(dat_pc_wob, 
                      aes(x=variable2, y=score)) +
  geom_jitter(height=0, color=dat_pc_wob$color,
              size=1.5, alpha=0.9)+   
  geom_violin(alpha=0.3) +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0, 1.05)) +
  scale_x_discrete(labels=x_labels[1:4]) + 
  scale_fill_manual(values=c("#AAAAAA", "#444444")) +
  theme_classic() + labs(x="Age group", y="Proportion correct answers \nin the instructions comprehension") + 
  theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm")) +
  theme(text = element_text(size=14), legend.position="none") +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14))

ann_text<-data.frame(x=c(1, 2, 3, 4), y=rep(1.05, 4), label=rep("***", 4))

wob_by_group <- ggplot(dat_pc_wob, 
                       aes(x=variable2, y=score)) +
  geom_jitter(height=0, color=dat_pc_wob$color, size=1.5, alpha=0.9)+   
  geom_violin(alpha=0.3) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  scale_x_discrete(labels=x_labels[1:4]) + 
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0, 1.05)) + 
  scale_fill_manual(values=c("#AAAAAA", "#444444")) +
  annotate("text", x=ann_text$x, y=ann_text$y, label=ann_text$label, size = 6) + 
  theme_classic() + labs(x="Age group", y="Proportion correct answers \nin the EV comparisons task") + 
  theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm")) +
  theme(text = element_text(size=14), legend.position="none") +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14))

grid.arrange(pc_by_group, wob_by_group, ncol=2)

## Fishing choices
# t-tests

# Overall fishing choices compared to chance and children vs adults
# based on what participants see
t.test(subset(dat, trial_number == 1 & age_group_coded == 0)$fishing_non_Z, 
       mu=0.5, var.equal=TRUE)
t.test(subset(dat, trial_number == 1 & age_group_coded == 1)$fishing_non_Z, 
       mu=0.5, var.equal=TRUE)
t.test(subset(dat, trial_number == 1 & age_group_coded == 2)$fishing_non_Z, 
       mu=0.5, var.equal=TRUE)
t.test(subset(dat, trial_number == 1 & age_group_coded == 3)$fishing_non_Z, 
       mu=0.5, var.equal=TRUE)

# normality tests
shapiro.test(subset(dat_4_5, trial_number ==1)$fishing_non_Z)
shapiro.test(subset(dat_6_7, trial_number ==1)$fishing_non_Z)
shapiro.test(subset(dat_8_9, trial_number ==1)$fishing_non_Z)
shapiro.test(subset(dat_10_12, trial_number ==1)$fishing_non_Z)

# non-parametric tests
wilcox.test(subset(dat_4_5, trial_number ==1)$fishing_non_Z, mu=0.5, exact=FALSE)
wilcox.test(subset(dat_6_7, trial_number ==1)$fishing_non_Z, mu=0.5, exact=FALSE)
wilcox.test(subset(dat_8_9, trial_number ==1)$fishing_non_Z, mu=0.5, exact=FALSE)
wilcox.test(subset(dat_10_12, trial_number ==1)$fishing_non_Z, mu=0.5, exact=FALSE)

subset(dat_4_5, trial_number ==1) %>% rstatix::wilcox_effsize(fishing_non_Z ~ 1, mu = 0.5)
subset(dat_6_7, trial_number ==1) %>% rstatix::wilcox_effsize(fishing_non_Z ~ 1, mu = 0.5)
subset(dat_8_9, trial_number ==1) %>% rstatix::wilcox_effsize(fishing_non_Z ~ 1, mu = 0.5)
subset(dat_10_12, trial_number ==1) %>% rstatix::wilcox_effsize(fishing_non_Z ~ 1, mu = 0.5)



# Plot average scores by age group
# Create dataset for fishing choices
age_group <- factor(c("4-5", "6-7", "8-9", "10-12", "18+"), levels=c("4-5", "6-7", "8-9", "10-12", "Adults"))
x_points <- c(0.1, 0.2, 0.3, 0.4, 0.5)
x_labels <- c("4-5", "6-7", "8-9", "10-12", "18+")

# based on what participants see
fishing <- c(
  mean(subset(dat, trial_number == 1 & age_group_coded == 0)$fishing_non_Z),
  mean(subset(dat, trial_number == 1 & age_group_coded == 1)$fishing_non_Z),
  mean(subset(dat, trial_number == 1 & age_group_coded == 2)$fishing_non_Z),
  mean(subset(dat, trial_number == 1 & age_group_coded == 3)$fishing_non_Z),
  mean(subset(adu_dat, trial_number == 1)$fishing_non_Z)
)
fishing_sem <- c(
  sem(subset(dat, trial_number == 1 & age_group_coded == 0)$fishing_non_Z),
  sem(subset(dat, trial_number == 1 & age_group_coded == 1)$fishing_non_Z),
  sem(subset(dat, trial_number == 1 & age_group_coded == 2)$fishing_non_Z),
  sem(subset(dat, trial_number == 1 & age_group_coded == 3)$fishing_non_Z),
  sem(subset(adu_dat, trial_number == 1)$fishing_non_Z)
)

fishing_choices_dat_by_group <- data.frame(x_points, x_labels,
                                           fishing, fishing_sem)

fishing_choices_by_group_plot <- ggplot(fishing_choices_dat_by_group) + 
  geom_line(aes(y=fishing, x=x_points), size=1, color="darkgrey")+
  geom_point(aes(y=fishing, x=x_points), color="darkgrey", size = 4)+
  geom_ribbon(aes(ymin=fishing-fishing_sem, 
                  ymax=fishing+fishing_sem, 
                  x=x_points), fill="darkgrey" , alpha = 0.3)+
  
  scale_x_continuous(labels=x_labels) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  theme_classic() +
  theme(text = element_text(size=14), legend.position="none") +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14)) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  labs(x="Age group", y="Correct fishing choices")

fishing_choices_by_group_plot_child <- ggplot(fishing_choices_dat_by_group[1:4,]) + 
  geom_line(aes(y=fishing, x=x_points), size=1, color="darkgrey")+
  geom_point(aes(y=fishing, x=x_points), color="darkgrey", size = 4)+
  geom_ribbon(aes(ymin=fishing-fishing_sem, 
                  ymax=fishing+fishing_sem, 
                  x=x_points), fill="darkgrey" , alpha = 0.3)+
  
  scale_x_continuous(labels=x_labels[1:4]) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)) +
  theme_classic() +
  theme(text = element_text(size=14), legend.position="none") +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14)) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  labs(x="Age group", y="Correct fishing choices")

## Violin plots

fishing <- c(
  (subset(dat, trial_number == 1 & age_group_coded == 0)$fishing_non_Z),
  (subset(dat, trial_number == 1 & age_group_coded == 1)$fishing_non_Z),
  (subset(dat, trial_number == 1 & age_group_coded == 2)$fishing_non_Z),
  (subset(dat, trial_number == 1 & age_group_coded == 3)$fishing_non_Z),
  (subset(adu_dat, trial_number == 1)$fishing_non_Z)
)
fishing_sem <- c(
  (subset(dat, trial_number == 1 & age_group_coded == 0)$fishing_non_Z),
  (subset(dat, trial_number == 1 & age_group_coded == 1)$fishing_non_Z),
  (subset(dat, trial_number == 1 & age_group_coded == 2)$fishing_non_Z),
  (subset(dat, trial_number == 1 & age_group_coded == 3)$fishing_non_Z),
  (subset(adu_dat, trial_number == 1)$fishing_non_Z)
)

x_points <- c(rep(0.1, length(subset(dat, trial_number == 1 & age_group_coded == 0)$fishing_non_Z)),
              rep(0.2, length(subset(dat, trial_number == 1 & age_group_coded == 1)$fishing_non_Z)),
              rep(0.3, length(subset(dat, trial_number == 1 & age_group_coded == 2)$fishing_non_Z)),
              rep(0.4, length(subset(dat, trial_number == 1 & age_group_coded == 3)$fishing_non_Z)),
              rep(0.5, length(subset(adu_dat, trial_number == 1)$fishing_non_Z)))

x_labels <- c(rep("4-5", length(subset(dat, trial_number == 1 & age_group_coded == 0)$fishing_non_Z)),
              rep("6-7", length(subset(dat, trial_number == 1 & age_group_coded == 1)$fishing_non_Z)),
              rep("8-9", length(subset(dat, trial_number == 1 & age_group_coded == 2)$fishing_non_Z)),
              rep("10-12", length(subset(dat, trial_number == 1 & age_group_coded == 3)$fishing_non_Z)),
              rep("18+", length(subset(adu_dat, trial_number == 1)$fishing_non_Z)))

group <- c(rep("Children", length(subset(dat, trial_number == 1)$fishing_non_Z)), 
           rep("Adults", length(subset(adu_dat, trial_number == 1)$fishing_non_Z)))

fishing_choices_dat_by_group <- data.frame(x_points, x_labels, group,
                                           fishing, fishing_sem)

ann_text<-data.frame(x=c(1, 2, 3, 4), y=rep(1.05, 4), label=rep("***", 4))

fishing_by_group <- ggplot(subset(fishing_choices_dat_by_group, group == "Children"), 
                             aes(x=as.factor(x_points), y=fishing, fill=x_labels)) +
  geom_jitter(height=0, color ="#BEBEBE", size=1.5, alpha=0.9)+   
  geom_violin(alpha=0.3) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  scale_x_discrete(labels=c("4-5", "6-7", "8-9", "10-12")) + 
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0, 1.05)) + 
  scale_fill_manual(values=c("#AAAAAA", "#AAAAAA", "#AAAAAA", "#AAAAAA")) +
  annotate("text", x=ann_text$x, y=ann_text$y, label=ann_text$label, size = 6) + 
  theme_classic() + labs(x="Age group", y="Proportion correct fishing choices") + 
  theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm")) +
  theme(text = element_text(size=14), legend.position="none") +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14))


#grid.arrange(pc_by_group, wob_by_group, fishing_choices_by_group_plot, ncol=3)
control_vars_plot <- grid.arrange(pc_by_group, wob_by_group, fishing_by_group, ncol=3)

#### Reaction times ####
dat$log_RT_info_choice_non_Z <- log(dat$RT_info_choice_non_Z)
mean_RT_info_choice_log <- stats::aggregate(dat$log_RT_info_choice_non_Z, list(dat$subject_ID), FUN=mean)$x
mean_RT_info_choice <- stats::aggregate(dat$RT_info_choice_non_Z, list(dat$subject_ID), FUN=mean)$x
mean_age_in_years_non_Z <- stats::aggregate(dat$age_in_years_non_Z, list(dat$subject_ID), FUN=mean)$x 
dat_rt <- data.frame(mean_RT_info_choice_log, mean_RT_info_choice, mean_age_in_years_non_Z)

cor.test(mean_RT_info_choice_log, mean_age_in_years_non_Z, method="pearson")
cor.test(mean_RT_info_choice, mean_age_in_years_non_Z, method="pearson")


#### Average fishing reward ####
##### With actual fishing choice #####
avg_reward_external <- aggregate(dat$reward_external, FUN=mean, by=list(dat$subject_ID))$x[order(unique(dat$subject_ID))]
age_in_years <- aggregate(dat$age_in_years_non_Z, FUN=mean, by=list(dat$subject_ID))$x[order(unique(dat$subject_ID))]

dat_rew <- data.frame(unique(dat$subject_ID),avg_reward_external, age_in_years)
cor.test(avg_reward_external, age_in_years, method="pearson")

#### Correlations between deltas ####
##### Matrices #####
corr_variables <- c("delta_EV", "delta_SD", "delta_uncertainty_level", "delta_agency")
dat_corr <- dat[corr_variables]
dat_corr_mat <- round(cor(dat_corr, method="pearson"), 3)
dat_p_mat <- ggcorrplot::cor_pmat(dat_corr)
ggcorrplot(dat_corr_mat, outline.col = "white", lab = TRUE, insig = "blank")
dat_corr_mat 
dat_p_mat

##### Plot example deltas over trial number #####
# Children
ex_child <- 111  # pick a child at random
ex_dat <- subset(dat, dat$subject_ID == (as.character(dat$subject_ID[ex_child])))
ex_delta_EV <- scale(ex_dat$delta_EV_non_Z, center=center_variables)
ex_delta_uncertainty_level <- scale(ex_dat$delta_uncertainty_level_non_Z, center=center_variables)
ex_delta_agency <- scale(ex_dat$delta_agency_non_Z, center=center_variables)
ex_deltas = data.frame(ex_delta_EV, ex_delta_uncertainty_level, ex_delta_agency)
ex_deltas$trial_number <- seq(1, length(ex_delta_agency))
ex_deltas <- tidyr::gather(ex_deltas, key="variable", value="score", ex_delta_EV:ex_delta_agency, factor_key=TRUE)

delta_values_plot <- ggplot(ex_deltas, aes(x=trial_number, y=score, color=variable))+
  geom_line(size=1) + scale_color_manual(values=c(ev_color, agn_color, unc_color), labels=c("ΔEV", "Δagency", "Δuncertainty level")) + 
  labs(y="Rescaled value", x="Trial number") +
  theme_classic2() +
  theme(text = element_text(size=16))
delta_values_plot


#### Models of children only ####
##### Models #####
# Age as a continuous variable interacting with other factors
{
  child_mod_full <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                       percent_comprehension + wob + fishing +
                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

  child_mod_null <- glmer(info_choice ~ (1 | subject_ID),
                          data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  child_mod_drop_Af_int <- glmer(info_choice ~ (delta_uncertainty_level + delta_agency)*age_in_years + delta_EV + 
                                   percent_comprehension + wob + fishing +
                                   (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                 data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_drop_C_int <- glmer(info_choice ~ (delta_EV + delta_agency)*age_in_years + delta_uncertainty_level + 
                                   percent_comprehension + wob + fishing +
                                   (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                 data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_drop_Ac_int <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level)*age_in_years + delta_agency + 
                                   percent_comprehension + wob + fishing +
                                   (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                 data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_drop_Af_fixed <- glmer(info_choice ~ (delta_uncertainty_level + delta_agency)*age_in_years + delta_EV:age_in_years +
                       percent_comprehension + wob + fishing +
                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_drop_C_fixed <- glmer(info_choice ~ (delta_EV + delta_agency)*age_in_years + delta_uncertainty_level:age_in_years +
                                     percent_comprehension + wob + fishing +
                                     (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                   data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_drop_Ac_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level)*age_in_years + delta_agency:age_in_years +
                                    percent_comprehension + wob + fishing +
                                    (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                  data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_drop_age_in_years_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency):age_in_years +
                       delta_EV + delta_uncertainty_level + delta_agency +
                       percent_comprehension + wob + fishing +
                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  child_mod_drop_pc_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years + 
                                    wob + fishing +
                                   (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                 data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_drop_wob_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years + 
                                   percent_comprehension + fishing +
                                   (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                 data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_drop_fishing_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years + 
                                   percent_comprehension + wob +
                                   (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                 data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
}

##### Significance tests #####
sjPlot::tab_model(child_mod_full, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
test_child_EV_int <- anova(child_mod_full, child_mod_drop_Af_int, test="Chisq")
test_child_uncertainty_int <- anova(child_mod_full, child_mod_drop_C_int, test="Chisq")
test_child_agency_int <- anova(child_mod_full, child_mod_drop_Ac_int, test="Chisq")
test_child_EV_fixed <- anova(child_mod_full, child_mod_drop_Af_fixed, test="Chisq")
test_child_uncertainty_fixed <- anova(child_mod_full, child_mod_drop_C_fixed, test="Chisq")
test_child_agency_fixed <- anova(child_mod_full, child_mod_drop_Ac_fixed, test="Chisq")
test_child_age_in_years_fixed <- anova(child_mod_full, child_mod_drop_age_in_years_fixed, test="Chisq")
if (use_covariates) {
  test_child_pc_fixed <- anova(child_mod_full, child_mod_drop_pc_fixed, test="Chisq")
  test_child_wob_fixed <- anova(child_mod_full, child_mod_drop_wob_fixed, test="Chisq")
  test_child_fishing_fixed <- anova(child_mod_full, child_mod_drop_fishing_fixed, test="Chisq")

}


##### Control models ##### 
###### Additional interactions #####
{
  child_mod_full_int3 <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                             delta_EV:delta_uncertainty_level:age_in_years + 
                             delta_EV:delta_agency:age_in_years +
                             delta_agency:delta_uncertainty_level:age_in_years +
                             percent_comprehension + wob + fishing +
                             (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                           data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  sjPlot::tab_model(child_mod_full_int3, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
  
  child_mod_int3 <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                 delta_EV:delta_uncertainty_level:age_in_years + 
                                 delta_EV:delta_agency:age_in_years +
                                 delta_agency:delta_uncertainty_level:age_in_years +
                                 percent_comprehension + wob + fishing +
                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_Af_int <- glmer(info_choice ~ (delta_uncertainty_level + delta_agency)*age_in_years +
                                 delta_EV + 
                                 delta_EV:delta_uncertainty_level:age_in_years + 
                                 delta_EV:delta_agency:age_in_years +
                                 delta_agency:delta_uncertainty_level:age_in_years +
                                 percent_comprehension + wob + fishing +
                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_C_int <- glmer(info_choice ~ (delta_EV + delta_agency)*age_in_years +
                                 delta_uncertainty_level + 
                                 delta_EV:delta_uncertainty_level:age_in_years + 
                                 delta_EV:delta_agency:age_in_years +
                                 delta_agency:delta_uncertainty_level:age_in_years +
                                 percent_comprehension + wob + fishing +
                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_Ac_int <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level)*age_in_years +
                                 delta_agency +
                                 delta_EV:delta_uncertainty_level:age_in_years + 
                                 delta_EV:delta_agency:age_in_years +
                                 delta_agency:delta_uncertainty_level:age_in_years +
                                 percent_comprehension + wob + fishing +
                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_Af_fixed <- glmer(info_choice ~ (delta_uncertainty_level + delta_agency)*age_in_years +
                                delta_EV:age_in_years +
                                 delta_EV:delta_uncertainty_level:age_in_years + 
                                 delta_EV:delta_agency:age_in_years +
                                 delta_agency:delta_uncertainty_level:age_in_years +
                                 percent_comprehension + wob + fishing +
                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_C_fixed <- glmer(info_choice ~ (delta_EV + delta_agency)*age_in_years +
                                              delta_uncertainty_level:age_in_years +
                                             delta_EV:delta_uncertainty_level:age_in_years + 
                                               delta_EV:delta_agency:age_in_years +
                                               delta_agency:delta_uncertainty_level:age_in_years +
                                               percent_comprehension + wob + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_Ac_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level)*age_in_years +
                                              delta_agency:age_in_years +
                                              delta_EV:delta_uncertainty_level:age_in_years + 
                                              delta_EV:delta_agency:age_in_years +
                                              delta_agency:delta_uncertainty_level:age_in_years +
                                              percent_comprehension + wob + fishing +
                                              (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                            data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_AfCage_int <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                 delta_EV:delta_agency:age_in_years +
                                 delta_agency:delta_uncertainty_level:age_in_years +
                                 percent_comprehension + wob + fishing +
                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_AfAcage_int <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                 delta_EV:delta_uncertainty_level:age_in_years + 
                                 delta_agency:delta_uncertainty_level:age_in_years +
                                 percent_comprehension + wob + fishing +
                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_CAcage_int <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                 delta_EV:delta_uncertainty_level:age_in_years + 
                                 delta_EV:delta_agency:age_in_years +
                                 percent_comprehension + wob + fishing +
                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_age_in_years_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                (delta_EV + delta_uncertainty_level + delta_agency):age_in_years +
                                 delta_EV:delta_uncertainty_level:age_in_years + 
                                 delta_EV:delta_agency:age_in_years +
                                 delta_agency:delta_uncertainty_level:age_in_years +
                                 percent_comprehension + wob + fishing +
                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  child_mod_int3_drop_pc_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                                 delta_EV:delta_uncertainty_level:age_in_years + 
                                                 delta_EV:delta_agency:age_in_years +
                                                 delta_agency:delta_uncertainty_level:age_in_years +
                                                 wob + fishing +
                                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_int3_drop_wob_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                          delta_EV:delta_uncertainty_level:age_in_years + 
                                          delta_EV:delta_agency:age_in_years +
                                          delta_agency:delta_uncertainty_level:age_in_years +
                                          percent_comprehension + fishing +
                                          (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                        data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  child_mod_int3_drop_fishing_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                          delta_EV:delta_uncertainty_level:age_in_years + 
                                          delta_EV:delta_agency:age_in_years +
                                          delta_agency:delta_uncertainty_level:age_in_years +
                                          percent_comprehension + wob +
                                          (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                        data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
}

test_child_int3_EV_int <- anova(child_mod_full_int3, child_mod_int3_drop_Af_int, test="Chisq")
test_child_int3_uncertainty_int <- anova(child_mod_full_int3, child_mod_int3_drop_C_int, test="Chisq")
test_child_int3_agency_int <- anova(child_mod_full_int3, child_mod_int3_drop_Ac_int, test="Chisq")
test_child_int3_EV_fixed <- anova(child_mod_full_int3, child_mod_int3_drop_Af_fixed, test="Chisq")
test_child_int3_uncertainty_fixed <- anova(child_mod_full_int3, child_mod_int3_drop_C_fixed, test="Chisq")
test_child_int3_agency_fixed <- anova(child_mod_full_int3, child_mod_int3_drop_Ac_fixed, test="Chisq")
test_child_int3_age_in_years_fixed <- anova(child_mod_full_int3, child_mod_int3_drop_age_in_years_fixed, test="Chisq")
test_child_int3_pc_fixed <- anova(child_mod_full_int3, child_mod_int3_drop_pc_fixed, test="Chisq")
test_child_int3_wob_fixed <- anova(child_mod_full_int3, child_mod_int3_drop_wob_fixed, test="Chisq")
test_child_int3_fishing_fixed <- anova(child_mod_full_int3, child_mod_int3_drop_fishing_fixed, test="Chisq")
test_child_int3_EV_uncertainty_age_int <- anova(child_mod_full_int3, child_mod_int3_drop_AfCage_int, test="Chisq")
test_child_int3_EV_agency_age_int <- anova(child_mod_full_int3, child_mod_int3_drop_AfAcage_int, test="Chisq")
test_child_int3_uncertainty_agency_age_int <- anova(child_mod_full_int3, child_mod_int3_drop_CAcage_int, test="Chisq")


###### Delta agency is 0 ######
## Effect of EV and uncertainty when agency is 0/not 0
# Age in months is log, interactions with age, delta agency == 0
{
  child_mod_delta_agency_0 <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level)*age_in_years +
                             percent_comprehension + wob + fishing +
                             (delta_EV + delta_uncertainty_level | subject_ID),
                           data = subset(dat, dat$delta_agency==0), family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  sjPlot::tab_model(child_mod_delta_agency_0, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
  
  
  child_mod_delta_agency_0_drop_Af_int <- glmer(info_choice ~ delta_uncertainty_level*age_in_years +
                                              delta_EV + percent_comprehension + wob + fishing +
                                              (delta_EV + delta_uncertainty_level | subject_ID),
                                            data = subset(dat, dat$delta_agency==0), family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_delta_agency_0_drop_C_int <- glmer(info_choice ~ delta_EV *age_in_years +
                                              delta_uncertainty_level + percent_comprehension + wob + fishing +
                                              (delta_EV + delta_uncertainty_level | subject_ID),
                                            data = subset(dat, dat$delta_agency==0), family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_delta_agency_0_drop_Af_fixed <- glmer(info_choice ~ (delta_uncertainty_level)*age_in_years +
                                                   delta_EV:age_in_years + percent_comprehension + wob + fishing +
                                                   (delta_EV + delta_uncertainty_level | subject_ID),
                                                 data = subset(dat, dat$delta_agency==0), family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_delta_agency_0_drop_C_fixed <- glmer(info_choice ~ (delta_EV)*age_in_years +
                                                   delta_uncertainty_level:age_in_years + percent_comprehension + wob + fishing +
                                                   (delta_EV + delta_uncertainty_level | subject_ID),
                                                 data = subset(dat, dat$delta_agency==0), family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_delta_agency_0_drop_age_in_years_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level):age_in_years +
                                       delta_EV + delta_uncertainty_level + percent_comprehension + wob + fishing +
                                      (delta_EV + delta_uncertainty_level | subject_ID),
                                    data = subset(dat, dat$delta_agency==0), family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  child_mod_delta_agency_0_drop_pc_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level)*age_in_years +
                                                   wob + fishing +
                                                   (delta_EV + delta_uncertainty_level | subject_ID),
                                                 data = subset(dat, dat$delta_agency==0), family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_delta_agency_0_drop_wob_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level)*age_in_years +
                                                   percent_comprehension + fishing +
                                                   (delta_EV + delta_uncertainty_level | subject_ID),
                                                 data = subset(dat, dat$delta_agency==0), family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_delta_agency_0_drop_fishing_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level)*age_in_years +
                                                   percent_comprehension + wob +
                                                   (delta_EV + delta_uncertainty_level | subject_ID),
                                                 data = subset(dat, dat$delta_agency==0), family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
}  

test_child_delta_ag0_EV_int <- anova(child_mod_delta_agency_0, child_mod_delta_agency_0_drop_Af_int, test="Chisq")
test_child_delta_ag0_uncertainty_int <- anova(child_mod_delta_agency_0, child_mod_delta_agency_0_drop_C_int, test="Chisq")
test_child_delta_ag0_EV_fixed <- anova(child_mod_delta_agency_0, child_mod_delta_agency_0_drop_Af_fixed, test="Chisq")
test_child_delta_ag0_uncertainty_fixed <- anova(child_mod_delta_agency_0, child_mod_delta_agency_0_drop_C_fixed, test="Chisq")
test_child_delta_ag0_age_in_years_fixed <- anova(child_mod_delta_agency_0, child_mod_delta_agency_0_drop_age_in_years_fixed, test="Chisq")
test_child_delta_ag0_pc_fixed <- anova(child_mod_delta_agency_0, child_mod_delta_agency_0_drop_pc_fixed, test="Chisq")
test_child_delta_ag0_wob_fixed <- anova(child_mod_delta_agency_0, child_mod_delta_agency_0_drop_wob_fixed, test="Chisq")
test_child_delta_ag0_fishing_fixed <- anova(child_mod_delta_agency_0, child_mod_delta_agency_0_drop_fishing_fixed, test="Chisq")


##### Previous choices #####
## Effect of previous choices
# Age in months is log, with previous choices interacting with age
{
  child_mod_prev_choice <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency + prev_info_choice)*age_in_years +
                                          percent_comprehension + wob + fishing +
                                          (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                        data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  sjPlot::tab_model(child_mod_prev_choice, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
  
  
  child_mod_prev_choice_drop_Af_int <- glmer(info_choice ~ (delta_uncertainty_level + delta_agency + prev_info_choice)*age_in_years +
                                                           delta_EV + percent_comprehension + wob + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_prev_choice_drop_C_int <- glmer(info_choice ~ (delta_EV + delta_agency + prev_info_choice)*age_in_years +
                                                          delta_uncertainty_level + percent_comprehension + wob + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_prev_choice_drop_Ac_int <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + prev_info_choice)*age_in_years +
                                                           delta_agency + percent_comprehension + wob + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_prev_choice_drop_prev_int <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                                             prev_info_choice + percent_comprehension + wob + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_prev_choice_drop_Af_fixed <- glmer(info_choice ~ (delta_uncertainty_level + delta_agency + prev_info_choice)*age_in_years +
                                                        delta_EV:age_in_years +
                                                        percent_comprehension + wob + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_prev_choice_drop_C_fixed <- glmer(info_choice ~ (delta_EV + delta_agency + prev_info_choice)*age_in_years +
                                               delta_uncertainty_level:age_in_years +
                                               percent_comprehension + wob + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_prev_choice_drop_Ac_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + prev_info_choice)*age_in_years +
                                               delta_agency:age_in_years +
                                               percent_comprehension + wob + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_prev_choice_drop_prev_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                                               prev_info_choice:age_in_years + 
                                                               percent_comprehension + wob + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_prev_choice_drop_age_in_years_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency + prev_info_choice):age_in_years +
                                                                delta_EV + delta_uncertainty_level + delta_agency + prev_info_choice +
                                                                percent_comprehension + wob + fishing +
                                                 (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  
  child_mod_prev_choice_drop_pc_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency + prev_info_choice)*age_in_years +
                                                wob + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_prev_choice_drop_wob_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency + prev_info_choice)*age_in_years +
                                               percent_comprehension + fishing +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_prev_choice_drop_fishing_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency + prev_info_choice)*age_in_years +
                                               percent_comprehension + wob +
                                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                             data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
}

test_child_prev_EV_int <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_Af_int, test="Chisq")
test_child_prev_uncertainty_int <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_C_int, test="Chisq")
test_child_prev_agency_int <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_Ac_int, test="Chisq")
test_child_prev_prev_int <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_prev_int, test="Chisq")
test_child_prev_EV_fixed <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_Af_fixed, test="Chisq")
test_child_prev_uncertainty_fixed <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_C_fixed, test="Chisq")
test_child_prev_agency_fixed <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_Ac_fixed, test="Chisq")
test_child_prev_prev_fixed <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_prev_fixed, test="Chisq")
test_child_prev_age_in_years_fixed <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_age_in_years_fixed, test="Chisq")
test_child_prev_pc_fixed <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_pc_fixed, test="Chisq")
test_child_prev_wob_fixed <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_wob_fixed, test="Chisq")
test_child_prev_fishing_fixed <- anova(child_mod_prev_choice, child_mod_prev_choice_drop_fishing_fixed, test="Chisq")

sjPlot::tab_model(child_mod_prev_choice, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)


##### Interaction plots #####
child_mod_int <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*scale(age_in_years_non_Z, center=center_age, scale=rescale_age) +
                                 percent_comprehension + wob + fishing +
                                 (delta_EV + delta_uncertainty_level + delta_agency
                                  | subject_ID),
                               data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

if (data_batch == 1) {
  modx_vals <- c(12, 11, 10, 9, 8, 7, 6, 5, 4)
} else {
  modx_vals <- c(11, 10, 9, 8, 7, 6, 5, 4)
}

int_plot_EV <- interact_plot(child_mod_int, pred = delta_EV, modx = age_in_years_non_Z, interval = TRUE, 
                             plot.points = FALSE, int.width=0.95, line.thickness=1, 
                             legend.main = "Age (years)", modx.values = modx_vals, 
                             colors="red", reverse.lty=FALSE, data=dat) + theme_classic2() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x="ΔEV (right-left)", y="P(choose right)") + 
  theme(legend.box.spacing = unit(0, "pt")) +
  theme(text = element_text(size=12)) +
  theme(legend.key.width= unit(60, "pt")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size=10)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))


int_plot_uncertainty <- interact_plot(child_mod_int, pred = delta_uncertainty_level, modx = age_in_years_non_Z, interval = TRUE, 
                                      plot.points = FALSE, int.width = 0.95, line.thickness = 1, 
                                      legend.main = "Age (years)",
                                      modx.values = modx_vals, colors="blue", data=dat) + theme_classic2() +
  labs(x="Δuncertainty (right-left)", y="P(choose right)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.box.spacing = unit(0, "pt")) +
  theme(text = element_text(size=12)) +
  theme(legend.key.width= unit(60, "pt")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size=10)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

int_plot_agency <- interact_plot(child_mod_int, pred = delta_agency, modx = age_in_years_non_Z, interval = TRUE, 
                                 plot.points = FALSE, int.width = 0.95, line.thickness = 1, 
                                 legend.main = "Age (years)",
                                 modx.values = modx_vals, colors="green", data=dat) + theme_classic2() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x="Δagency (right-left)", y="P(choose right)") + 
  theme(legend.box.spacing = unit(0, "pt")) +
  theme(text = element_text(size=12)) +
  theme(legend.key.width= unit(60, "pt")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size=10)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

##### Johnson-Neyman plots #####
# estimate predictive effects
child_mod_int2 <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                         percent_comprehension + wob + fishing +
                         (delta_EV + delta_uncertainty_level + delta_agency
                          | subject_ID),
                       data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)


sig_color <- "black"
nonsig_color <- "grey"

n_tests <- (12*max(dat$age_in_years_non_Z)) - (12*min(dat$age_in_years_non_Z))
slopes_EV <- modelbased::estimate_slopes(child_mod_int2, trend = "delta_EV", at = "age_in_years", length = n_tests)
slopes_uncertainty_level <- modelbased::estimate_slopes(child_mod_int2, trend = "delta_uncertainty_level", at = "age_in_years", length = n_tests)
slopes_agency <- modelbased::estimate_slopes(child_mod_int2, trend = "delta_agency", at = "age_in_years", length = n_tests)

x_pos <- c(4, 5, 6, 7, 8, 9, 10, 11, 12)
x_pos_Z <- find_Z(x_pos, dat$age_in_years_non_Z)
lims <- find_Z(c(4,12), dat$age_in_years_non_Z)

# plot significance of the effect across time (Johnson-Neyman plot)
# modelbased package
jn_plot_EV <- plot(slopes_EV) + 
  scale_fill_manual(values=c(ev_color)) + 
  theme_classic2() +
  scale_y_continuous(limits = c(-0.5, 2)) +
  scale_x_continuous(breaks=x_pos_Z, labels=x_pos, limits=lims) +
  labs(x="Age (years)", y="Estimated *β* relating ΔEV<br>to information-seeking choices", title="") + 
  theme(text = element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.y = ggtext::element_markdown())

jn_plot_uncertainty <- plot(slopes_uncertainty_level) + 
  scale_fill_manual(values=c(nonsig_color, unc_color)) + 
  theme_classic2() +
  scale_y_continuous(limits = c(-0.5, 2)) +
  scale_x_continuous(breaks=x_pos_Z, labels=x_pos, limits=lims) +
  labs(x="Age (years)", y="Estimated *β* relating Δuncertainty<br>to information-seeking choices", title="") + 
  theme(text = element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.y = ggtext::element_markdown())

jn_plot_agency <- plot(slopes_agency) + 
  scale_fill_manual(values=c(agn_color)) +  
  theme_classic2() +
  scale_y_continuous(limits = c(-0.5, 2)) +
  scale_x_continuous(breaks=x_pos_Z, labels=x_pos, limits=lims) +
  labs(x="Age (years)", y="Estimated *β* relating Δagency<br>to information-seeking choices", title="") + 
  theme(text = element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.y = ggtext::element_markdown())

grid.arrange(grobs=list(jn_plot_EV, jn_plot_agency,jn_plot_uncertainty))



#### Individual betas #### 
id <- c()
for (i in unique(dat$subject_ID)) {id <- append(id, i)}
age_group <- c()
age_in_years <- c()
pc <- c()
wob <- c()
fishing <- c()
coeff_intercept_sbj <- c()
coeff_EV_sbj<-c()
coeff_uncertainty_sbj <- c()
coeff_agency_sbj <- c()

exclude_ids <- c()
center_sbj_deltas <- FALSE

for (i in seq_along(unique(dat$subject_ID))) {
  id_sbj <- unique(dat$subject_ID)[i]
  # organize in a data frame for subject = id_sbj
  dat_s <- subset(dat, dat$subject_ID == id_sbj)

  # update model-independent values
  age_group[i] <- dat_s$age_group[1]
  age_in_years[i] <- dat_s$age_in_years_non_Z[1]
  pc[i] <- dat_s$percent_comprehension_non_Z[1]
  wob[i] <- dat_s$wob_non_Z[1]
  fishing[i] <- dat_s$fishing_non_Z[1]

  # logistic regression
  dat_s$delta_EV <- scale(dat_s$delta_EV_non_Z, center=center_sbj_deltas)
  dat_s$delta_uncertainty_level <- scale(dat_s$delta_uncertainty_level_non_Z, center=center_sbj_deltas)
  dat_s$delta_agency <- scale(dat_s$delta_agency_non_Z, center=center_sbj_deltas)
  
  x <- as.matrix(dat_s[, c('delta_EV', 'delta_uncertainty_level', 'delta_agency')])
  y <- as.matrix(dat_s$info_choice)
  
  tryCatch(
    {
      # alpha = 0: ridge regression
      model <- glmnet(x=x, y=y, family = "binomial", alpha=0)
      coefficients <- coef(model, FALSE)
      coeff_intercept_sbj[i] <- coefficients[1]
      coeff_EV_sbj[i] <- coefficients[2]
      coeff_uncertainty_sbj[i] <- coefficients[3]
      coeff_agency_sbj[i] <- coefficients[4]
    },
    error=function(cond){
      # if the model does not converge (all ys are the same)
      coeff_intercept_sbj[i] <- NA
      coeff_EV_sbj[i] <- NA
      coeff_uncertainty_sbj[i] <- NA
      coeff_agency_sbj[i] <- NA
      message(paste("Error at participant", as.character(id_sbj)))
      message(cond)
      message("\nCoefficients will be set at NA")
      }
  )
  
}

# create dataframe
dat_coeffs <- data.frame(id, age_group, age_in_years, 
                         pc, wob, fishing,
                         coeff_intercept_sbj, coeff_EV_sbj, coeff_uncertainty_sbj, coeff_agency_sbj)


##### Correlations #####
# EV and age
cor.test(dat_coeffs$coeff_EV_sbj, dat_coeffs$age_in_years, method="pearson")
correlationBF(dat_coeffs$coeff_EV_sbj, dat_coeffs$age_in_years)

# Uncertainty and age
cor.test(dat_coeffs$coeff_uncertainty_sbj,dat_coeffs$age_in_years, method="pearson")

# Agency and age
cor.test(dat_coeffs$coeff_agency_sbj, dat_coeffs$age_in_years, method="pearson")

# Agency and fishing choices
cor.test(dat_coeffs$coeff_agency_sbj, dat_coeffs$fishing, method="pearson")


#### Individual betas on trials with 0 delta agency #### 
dat_delta_ag_0 <- subset(dat, dat$delta_agency==0)
id <- c()
for (i in unique(dat$subject_ID)) {id <- append(id, i)}
age_group <- c()
age_in_years <- c()
pc <- c()
wob <- c()
fishing <- c()
coeff_intercept_sbj <- c()
coeff_EV_sbj<-c()
coeff_uncertainty_sbj <- c()
exclude_ids <- c()
center_sbj_deltas <- FALSE

for (i in seq_along(unique(dat$subject_ID))) {
  id_sbj <- unique(dat$subject_ID)[i]
  # organize in a data frame for subject = id_sbj
  dat_s <- subset(dat_delta_ag_0, dat_delta_ag_0$subject_ID == id_sbj)
  
  # update model-independent values
  age_group[i] <- dat_s$age_group[1]
  age_in_years[i] <- dat_s$age_in_years_non_Z[1]
  pc[i] <- dat_s$percent_comprehension_non_Z[1]
  wob[i] <- dat_s$wob_non_Z[1]
  fishing[i] <- dat_s$fishing_non_Z[1]
  
  
  # logistic regression
  dat_s$delta_EV <- scale(dat_s$delta_EV_non_Z, center=center_sbj_deltas)
  dat_s$delta_uncertainty_level <- scale(dat_s$delta_uncertainty_level_non_Z, center=center_sbj_deltas)

  x <- as.matrix(dat_s[, c('delta_EV', 'delta_uncertainty_level')])
  y <- as.matrix(dat_s$info_choice)
  
  tryCatch(
    {
      # alpha = 0: ridge regression
      model <- glmnet(x=x, y=y, family = "binomial", alpha=0)
      coefficients <- coef(model, FALSE)
      coeff_intercept_sbj[i] <- coefficients[1]
      coeff_EV_sbj[i] <- coefficients[2]
      coeff_uncertainty_sbj[i] <- coefficients[3]
    },
    error=function(cond){
      # if the model does not converge (all ys are the same)
      coeff_intercept_sbj[i] <- NA
      coeff_EV_sbj[i] <- NA
      coeff_uncertainty_sbj[i] <- NA
      message(paste("Error at participant", as.character(id_sbj)))
      message(cond)
      message("\nCoefficients will be set at NA")
    }
  )
  
}

# create dataframe
dat_coeffs_delta_ag_0 <- data.frame(id, age_group, age_in_years, 
                         pc, wob, fishing,
                         coeff_intercept_sbj, coeff_EV_sbj, coeff_uncertainty_sbj)


# EV and age
cor.test(dat_coeffs_delta_ag_0$coeff_EV_sbj, dat_coeffs_delta_ag_0$age_in_years, method="pearson")
correlationBF(dat_coeffs_delta_ag_0$coeff_EV_sbj, dat_coeffs_delta_ag_0$age_in_years)

# Uncertainty and age
cor.test(dat_coeffs_delta_ag_0$coeff_uncertainty_sbj,dat_coeffs_delta_ag_0$age_in_years, method="pearson")
