# AUTHOR: GAIA MOLINARO
# Data analysis for the main and supplementary experiments
# in Molinaro, Cogliati Dezza, & Sharot (in prep.)

# SET WORKING DIRECTORY
# wd = "insert_your_wd_here"
wd = "C:/Gaia/ABL"

#### LOAD REQUIRED PACKAGES #####
packages <- c("ggplot2", "GGally", "ggpubr", "afex", "plyr", "psych", "nls2", 
              "grid", "simr", "optimx", "plotrix", "BayesFactor", "ggcorrplot", 
              "tidyverse", "glue", "sjPlot", "sjmisc", "gridExtra", "rsq", 
              "sets", "qpcR")
pacman::p_load(packages, character.only = TRUE)


#### MAIN EXPERIMENT (COMPETITION) ####
#### Load and prep datasets #### 
setwd("C:/Gaia/ABL")
dat_compt <- read.csv("lookit_competition_check5_trials.csv")
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
dat_compt$non_z_delta_EV <- (as.numeric(dat_compt$delta_EV))
dat_compt$non_z_delta_uncertainty_level <- (as.numeric(dat_compt$delta_uncertainty_level))
dat_compt$non_z_delta_agency<- (as.numeric(dat_compt$delta_agency))
dat_compt$delta_EV <- scale(as.numeric(dat_compt$delta_EV))
dat_compt$delta_uncertainty_level <- scale(as.numeric(dat_compt$delta_uncertainty_level))
dat_compt$delta_IU <- scale(as.numeric(dat_compt$delta_IU))
dat_compt$delta_agency <- scale(as.numeric(dat_compt$delta_agency))
dat_compt$age_in_years <- scale(as.numeric(dat_compt$age_in_years))
dat_compt$age_group_coded <- as.factor(dat_compt$age_group_coded)
dat_compt$RT_info_choice <- scale(as.numeric(dat_compt$RT_info_choice))
dat_compt$gorilla_ID <- as.factor(dat_compt$gorilla_ID)
dat_compt$gender_coded <- as.factor(dat_compt$gender_coded)
contrasts(dat_compt$gender_coded) <- contr.helmert(length(table(dat_compt$gender_coded)))
dat_compt$info_choice <- as.factor(dat_compt$info_choice) # 0 = left, 1 = right
dat_compt$percent_comprehension_non_Z <- as.numeric(dat_compt$percent_comprehension)
dat_compt$percent_comprehension <- scale(as.numeric(dat_compt$percent_comprehension))
dat_compt$wob_non_Z <- as.numeric(dat_compt$wob)
dat_compt$wob <- scale(as.numeric(dat_compt$wob))
dat_compt$chance <- 0.5

# Set subject IDs
N <- length(sets::as.set(dat_compt$gorilla_ID))
tot_trials <- length(sets::as.set(dat_compt$trial_number))
subject_ID <- c()
for (n in 1:N) {
  for (t in 1:tot_trials){
    subject_ID <- c(subject_ID, glue("subject_{n}"))
  }
}
dat_compt$subject_ID <- as.factor(subject_ID)

# Create subsets of data for each age group
dat_compt_4_5 <- subset(dat_compt, age_group == "4-5")
dat_compt_6_7 <- subset(dat_compt, age_group == "6-7")
dat_compt_8_9 <- subset(dat_compt, age_group == "8-9")
dat_compt_10_12 <- subset(dat_compt, age_group == "10-12")

## Competition between motives - adults
setwd(wd)
adu_dat_compt <- read.csv("adults_competition_check1_trials.csv")
adu_dat_compt <- subset(adu_dat_compt, !(condition %in% c("catch_1", "catch_2")))
adu_dat_compt <- subset(adu_dat_compt, (catch_trials_score == 100) & (wob > 0.5))

# Fix trial number
adu_dat_compt$trial_number <- substr(adu_dat_compt$condition, 7, 8)

# Turn variables into appropriate data types
adu_dat_compt$delta_EV <- scale(as.numeric(adu_dat_compt$delta_EV))
adu_dat_compt$delta_uncertainty_level <- scale(as.numeric(adu_dat_compt$delta_uncertainty_level))
adu_dat_compt$delta_IU <- scale(as.numeric(adu_dat_compt$delta_IU))
adu_dat_compt$delta_agency <- scale(as.numeric(adu_dat_compt$delta_agency))
adu_dat_compt$age_in_years_non_z <- as.numeric(adu_dat_compt$age_in_years)
adu_dat_compt$age_in_years <- scale(as.numeric(adu_dat_compt$age_in_years))
adu_dat_compt$age_group <- "Adult"
adu_dat_compt$age_group_coded <- as.factor(4)
adu_dat_compt$RT_info_choice <- scale(as.numeric(adu_dat_compt$RT_info_choice))
adu_dat_compt$gorilla_ID <- as.factor(adu_dat_compt$gorilla_ID)
adu_dat_compt$gender_coded <- as.factor(adu_dat_compt$gender_coded)
contrasts(adu_dat_compt$gender_coded) <- contr.helmert(3)
adu_dat_compt$info_choice <- as.factor(adu_dat_compt$info_choice) # 0 = left, 1 = right
adu_dat_compt$percent_comprehension_non_Z <- as.numeric(adu_dat_compt$percent_comprehension)
adu_dat_compt$percent_comprehension <- scale(as.numeric(adu_dat_compt$percent_comprehension))
adu_dat_compt$wob_non_Z <- as.numeric(adu_dat_compt$wob)
adu_dat_compt$wob <- scale(as.numeric(adu_dat_compt$wob))
adu_dat_compt$chance <- 0.5

# Set subject IDs
adu_N <- length(sets::as.set(adu_dat_compt$gorilla_ID))
adu_tot_trials <- length(sets::as.set(adu_dat_compt$trial_number))
adu_subject_ID <- c()
for (n in 1:adu_N) {
  for (t in 1:adu_tot_trials){
    adu_subject_ID <- c(adu_subject_ID, glue("subject_{n}"))
  }
}
adu_dat_compt <- as.factor(adu_subject_ID)

## Merged datasets (children + adults)
# Create data subsets
cols <- c("gorilla_ID", "info_choice", 
          "delta_EV", "delta_uncertainty_level", "delta_agency", 
          "gender_coded", "wob", "percent_comprehension", 
          "wob_non_Z", "percent_comprehension_non_Z", 
          "age_group", "age_group_coded")
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

# Fix gender contrasts (children only had 2 levels, adults had 3)
contrasts(dat_compt_child_adu$gender_coded) <- contr.helmert(3)
contrasts(dat_compt_4_5_adu$gender_coded) <- contr.helmert(3)
contrasts(dat_compt_6_7_adu$gender_coded) <- contr.helmert(3)
contrasts(dat_compt_8_9_adu$gender_coded) <- contr.helmert(3)
contrasts(dat_compt_10_12_adu$gender_coded) <- contr.helmert(3)

# Create new variable and  contrasts for children vs adults factor
dat_compt_child_adu$group <- as.factor(ifelse(dat_compt_child_adu$age_group_coded == 4, "Adults", "Children"))
dat_compt_4_5_adu$group <- as.factor(ifelse(dat_compt_4_5_adu$age_group_coded == 4, "Adults", "Children"))
dat_compt_6_7_adu$group <-  as.factor(ifelse(dat_compt_6_7_adu$age_group_coded == 4, "Adults", "Children"))
dat_compt_8_9_adu$group <-  as.factor(ifelse(dat_compt_8_9_adu$age_group_coded == 4, "Adults", "Children"))
dat_compt_10_12_adu$group <-  as.factor(ifelse(dat_compt_10_12_adu$age_group_coded == 4, "Adults", "Children"))
contrasts(dat_compt_child_adu$group) <- contr.helmert(2)
contrasts(dat_compt_4_5_adu$group) <- contr.helmert(2)
contrasts(dat_compt_6_7_adu$group) <- contr.helmert(2)
contrasts(dat_compt_8_9_adu$group) <- contr.helmert(2)
contrasts(dat_compt_10_12_adu$group) <- contr.helmert(2)

#### Comprehension and wob scores ####
## t-tests

# Within group
t.test(subset(dat_compt, trial_number ==1)$wob_non_Z, mu=0.5)
t.test(subset(adu_dat_compt, trial_number == 1)$wob_non_Z, mu=0.5)

# Between groups
t.test(subset(dat_compt, trial_number ==1)$percent_comprehension_non_Z, 
       subset(adu_dat_compt, trial_number == 1)$percent_comprehension_non_Z)
t.test(subset(dat_compt, trial_number ==1)$wob_non_Z, 
       subset(adu_dat_compt, trial_number == 1)$wob_non_Z)

## Create long form data sets

# Create separate data sets with variables for children and adults' pc and wob
dat_compt_pc_wob <- dplyr::select(subset(dat_compt, trial_number==1), 
                             all_of(c("gorilla_ID", "age_group_coded", 
                                      "percent_comprehension_non_Z", 
                                      "wob_non_Z")))
dat_compt_pc_wob$percent_comprehension_child <- dat_compt_pc_wob$"percent_comprehension_non_Z"
dat_compt_pc_wob$wob_child <- dat_compt_pc_wob$"wob_non_Z"
adu_dat_compt_pc_wob <- dplyr::select(subset(adu_dat_compt, trial_number==1), 
                                  all_of(c("gorilla_ID", "age_group_coded", 
                                           "percent_comprehension_non_Z", 
                                           "wob_non_Z")))

# Merge datasets
dat_compt_child_adu_pc_wob <- rbind(dplyr::select(dat_compt_pc_wob, all_of(
  c("gorilla_ID", "age_group_coded", "percent_comprehension_non_Z", "wob_non_Z"))), 
  dplyr::select(adu_dat_compt_pc_wob, all_of(
    c("gorilla_ID", "age_group_coded", "percent_comprehension_non_Z", "wob_non_Z"))))

# Create group factor 
dat_compt_child_adu_pc_wob$group <- as.factor(ifelse(dat_compt_child_adu_pc_wob$age_group_coded == 4, "Adults", "Children"))

# Turn into long form 
dat_compt_child_adu_pc_wob <- dat_compt_child_adu_pc_wob %>%
  tidyr::gather(key = "variable", value = "score", percent_comprehension_non_Z:wob_non_Z, factor_key=TRUE)

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
  theme_classic() + labs(x="comprehension", y="Score") + theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm"))
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
  theme_classic() + labs(x="comprehension", y="Score") + theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm"))
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
  theme_classic() + labs(x="comprehension", y="Score") + theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm"))
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
  theme_classic() + labs(x="comprehension", y="Score") + theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm"))
compt_wob_by_group

grid.arrange(compt_pc_by_group, compt_wob_by_group, ncol=2)

####  Models by age group #### 
# 4-5
# Chance only
compt_mod_4_5_chance <- glmer(info_choice ~ 0 + chance + (1|gorilla_ID), 
                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full
compt_mod_4_5_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_uncertainty_level + delta_agency 
                               | gorilla_ID), 
                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect only
compt_mod_4_5_Af <- glmer(info_choice ~ delta_EV +
                              percent_comprehension + wob + gender_coded + 
                            (delta_EV 
                             | gorilla_ID), 
                          data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition only
compt_mod_4_5_C <- glmer(info_choice ~  delta_uncertainty_level +
                              percent_comprehension + wob + gender_coded + 
                            (delta_uncertainty_level 
                            | gorilla_ID), 
                         data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Action only
compt_mod_4_5_Ac <- glmer(info_choice ~ delta_agency +
                              percent_comprehension + wob + gender_coded +
                             (delta_agency 
                             | gorilla_ID), 
                          data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Cognition
compt_mod_4_5_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_uncertainty_level 
                              | gorilla_ID), 
                           data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Action 
compt_mod_4_5_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                               (delta_EV + delta_agency
                               | gorilla_ID), 
                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition + Action 
compt_mod_4_5_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                              (delta_uncertainty_level + delta_agency 
                              | gorilla_ID), 
                            data = dat_compt_4_5, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# 6-7
# Chance only
compt_mod_6_7_chance <- glmer(info_choice ~ 0 + chance + (1|gorilla_ID), 
                              data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
                             
# Full 
compt_mod_6_7_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_uncertainty_level + delta_agency                               | gorilla_ID), 
                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
                           
# Affect only
compt_mod_6_7_Af <- glmer(info_choice ~ delta_EV +
                            percent_comprehension + wob + gender_coded + 
                            (delta_EV
                             | gorilla_ID), 
                          data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition only 
compt_mod_6_7_C <- glmer(info_choice ~  delta_uncertainty_level +
                           percent_comprehension + wob + gender_coded + 
                           (delta_uncertainty_level 
                            | gorilla_ID), 
                        data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Action only
compt_mod_6_7_Ac <- glmer(info_choice ~ delta_agency +
                            percent_comprehension + wob + gender_coded +
                            (delta_agency 
                             | gorilla_ID), 
                          data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Cognition 
compt_mod_6_7_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + 
                             (delta_EV + delta_uncertainty_level 
                              | gorilla_ID), 
                           data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Action 
compt_mod_6_7_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_agency
                               | gorilla_ID), 
                            data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition + Action 
compt_mod_6_7_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                             (delta_uncertainty_level + delta_agency 
                              | gorilla_ID), 
                           data = dat_compt_6_7, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# 8-9
# Chance only
compt_mod_8_9_chance <- glmer(info_choice ~ 0 + chance + (1|gorilla_ID), 
                              data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
                                                           
# Full 
compt_mod_8_9_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_uncertainty_level + delta_agency
                               | gorilla_ID), 
                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
compt_mod_8_9_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                             (delta_EV + delta_uncertainty_level delta_agency +
                              | gorilla_ID), 
                           data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect only
compt_mod_8_9_Af <- glmer(info_choice ~ delta_EV +
                            percent_comprehension + wob + gender_coded + 
                            (delta_EV 
                             | gorilla_ID), 
                          data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition only
compt_mod_8_9_C <- glmer(info_choice ~  delta_uncertainty_level +
                           percent_comprehension + wob + gender_coded + 
                           (delta_uncertainty_level 
                            | gorilla_ID), 
                         data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Action only
compt_mod_8_9_Ac <- glmer(info_choice ~ delta_agency +
                            percent_comprehension + wob + gender_coded +
                            (delta_agency 
                             | gorilla_ID), 
                          data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Cognition
compt_mod_8_9_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + 
                             (delta_EV + delta_uncertainty_level 
                              | gorilla_ID), 
                          data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Action 
compt_mod_8_9_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_agency
                               | gorilla_ID), 
                            data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition + Action 
compt_mod_8_9_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                             (delta_uncertainty_level + delta_agency 
                              | gorilla_ID), 
                           data = dat_compt_8_9, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# 10-12
# Chance only
compt_mod_10_12_chance <- glmer(info_choice ~ 0 + chance + (1|gorilla_ID), 
                                data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full 
compt_mod_10_12_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_uncertainty_level + delta_agency +
                               | gorilla_ID), 
                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect only 
compt_mod_10_12_Af <- glmer(info_choice ~ delta_EV +
                            percent_comprehension + wob + gender_coded + 
                            (delta_EV 
                             | gorilla_ID), 
                            data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition only
compt_mod_10_12_C <- glmer(info_choice ~  delta_uncertainty_level +
                           percent_comprehension + wob + gender_coded + 
                           (delta_uncertainty_level 
                            | gorilla_ID), 
                          data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Action only 
compt_mod_10_12_Ac <- glmer(info_choice ~ delta_agency +
                            percent_comprehension + wob + gender_coded +
                            (delta_agency 
                             | gorilla_ID), 
                            data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Cognition 
compt_mod_10_12_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + 
                             (delta_EV + delta_uncertainty_level 
                              | gorilla_ID), 
                             data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Action 
compt_mod_10_12_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_agency
                               | gorilla_ID), 
                              data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition + Action 
compt_mod_10_12_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                             (delta_uncertainty_level + delta_agency 
                              | gorilla_ID), 
                            data = dat_compt_10_12, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Adults
# Chance only
adu_compt_mod_chance <- glmer(info_choice ~ 0 + chance + (1|gorilla_ID), 
                              data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full 
adu_compt_mod_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_uncertainty_level + delta_agency 
                               | gorilla_ID), 
                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect only
adu_compt_mod_Af <- glmer(info_choice ~ delta_EV +
                            percent_comprehension + wob + gender_coded + 
                            (delta_EV 
                             | gorilla_ID), 
                          data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition only
adu_compt_mod_C <- glmer(info_choice ~  delta_uncertainty_level +
                           percent_comprehension + wob + gender_coded + 
                           (delta_uncertainty_level 
                            | gorilla_ID), 
                         data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Action only 
adu_compt_mod_Ac <- glmer(info_choice ~ delta_agency +
                            percent_comprehension + wob + gender_coded +
                            (delta_agency 
                             | gorilla_ID), 
                          data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Cognition 
adu_compt_mod_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + 
                             (delta_EV + delta_uncertainty_level 
                              | gorilla_ID), 
                           data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Action 
adu_compt_mod_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_agency
                               | gorilla_ID), 
                            data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition + Action 
adu_compt_mod_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                             (delta_uncertainty_level + delta_agency 
                              | gorilla_ID), 
                           data = adu_dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
                         
# Children overall 
# Chance only
child_compt_mod_chance <- glmer(info_choice ~ 0 + chance + (1|gorilla_ID),
                                data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, with interactions
child_compt_mod_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + gender_coded + age_in_years +
                                  age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                  (delta_EV + delta_uncertainty_level + delta_agency
                                     | gorilla_ID), 
                                data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect only
child_compt_mod_Af <- glmer(info_choice ~ delta_EV +
                              percent_comprehension + wob + gender_coded + age_in_years +
                              delta_EV:age_in_years +
                              (delta_EV 
                               | gorilla_ID), 
                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition only 
child_compt_mod_C <- glmer(info_choice ~  delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + age_in_years +
                             delta_uncertainty_level:age_in_years +
                             (delta_uncertainty_level 
                              | gorilla_ID),
                           data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Action only
child_compt_mod_Ac <- glmer(info_choice ~ delta_agency +
                              percent_comprehension + wob + gender_coded + age_in_years +
                              delta_agency:age_in_years +
                              (delta_agency 
                               | gorilla_ID), 
                              data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Cognition 
child_compt_mod_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               percent_comprehension + wob + gender_coded + age_in_years +
                               delta_EV:age_in_years + delta_uncertainty_level:age_in_years +
                               (delta_EV + delta_uncertainty_level 
                                | gorilla_ID),
                             data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Action 
child_compt_mod_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                percent_comprehension + wob + gender_coded + age_in_years +
                                delta_EV:age_in_years + delta_agency:age_in_years +
                                (delta_EV + delta_agency
                                 | gorilla_ID),
                            data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition +  Action
child_compt_mod_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + age_in_years +
                               delta_uncertainty_level:age_in_years + delta_agency:age_in_years +
                               (delta_uncertainty_level + delta_agency 
                                | gorilla_ID),
                             data = dat_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
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
  theme(text = element_text(size=15), legend.position="none") +
  scale_x_continuous(labels=compt_betas_dat_by_group$x_labels) +
  scale_y_continuous(limits=c(-0.5, 2.25), breaks=seq(-0.5, 2.25, 0.25)) +
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
  scale_y_continuous(limits=c(-0.5, 2.25), breaks=seq(-0.5, 2.25, 0.25)) +
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
  scale_y_continuous(limits=c(-0.5, 2.25), breaks=seq(-0.5, 2.25, 0.25)) +
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
  scale_y_continuous(limits=c(-1.75, 2.25), breaks=seq(-1.75, 2.25, 0.25)) +
  labs(x="Age group", y="Standardized beta coefficient Intercept") +
  geom_hline(yintercept=0)
#+ geom_vline(xintercept=0.45, linetype="dotted")
compt_intercept_betas_by_group

grid.arrange(compt_hedonic_betas_by_group, compt_cognitive_betas_by_group, compt_instrumental_betas_by_group, ncol=3)

####  Children vs adults tests #### 
# Models
compt_mod_child_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + group + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | gorilla_ID), 
                           data = dat_compt_child_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
compt_mod_4_5_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
    percent_comprehension + wob + gender_coded + group + 
    group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
    (delta_EV + delta_uncertainty_level + delta_agency | gorilla_ID), 
    data = dat_compt_4_5_adu, family = binomial, 
    control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
compt_mod_6_7_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + group + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | gorilla_ID), 
                           data = dat_compt_6_7_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
compt_mod_8_9_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + group + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | gorilla_ID), 
                           data = dat_compt_8_9_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
compt_mod_10_12_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + group + 
                             group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                             (delta_EV + delta_uncertainty_level + delta_agency | gorilla_ID), 
                           data = dat_compt_10_12_adu, family = binomial, 
                           control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Model results
sjPlot::tab_model(compt_mod_child_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
sjPlot::tab_model(compt_mod_4_5_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
sjPlot::tab_model(compt_mod_6_7_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
sjPlot::tab_model(compt_mod_8_9_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE,show.ci=FALSE, show.se=TRUE)
sjPlot::tab_model(compt_mod_10_12_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)
