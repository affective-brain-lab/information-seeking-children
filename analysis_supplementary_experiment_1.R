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


#### SUPPLEMENTARY EXPERIMENT 1 (NO COMPETITION) ####
setwd(wd)
dat_no_compt <- read.csv("data_supplementary_experiment_1_children.csv")
dat_no_compt <- subset(dat_no_compt, !(condition %in% c("catch_1", "catch_2")))
dat_no_compt <- subset(dat_no_compt, !(age_in_years == 13))
dat_no_compt <- subset(dat_no_compt, catch_trials_score == 100)
dat_no_compt <- subset(dat_no_compt, wob > 0.5)
# Rename variables with mislabelled heading
dat_no_compt$subject_ID <- dat_no_compt$ï..subject_ID

# Create age groups
dat_no_compt$age_group <- "None"
dat_no_compt$age_group[dat_no_compt$age_in_years  %in% c(4, 5)] <- "4-5"
dat_no_compt$age_group[dat_no_compt$age_in_years  %in% c(6, 7)] <- "6-7"
dat_no_compt$age_group[dat_no_compt$age_in_years  %in% c(8, 9)] <- "8-9"
dat_no_compt$age_group[dat_no_compt$age_in_years  %in% c(10, 11, 12)] <- "10-12"
dat_no_compt$age_group <- factor(dat_no_compt$age_group, levels=c("4-5", "6-7", "8-9","10-12"))

# Turn variables into appropriate data types
# Use scale() for standardized values
dat_no_compt$delta_EV_non_Z <- (as.numeric(dat_no_compt$delta_EV))
dat_no_compt$delta_uncertainty_level_non_Z <- (as.numeric(dat_no_compt$delta_uncertainty_level))
dat_no_compt$delta_agency_non_Z<- (as.numeric(dat_no_compt$delta_agency))
dat_no_compt$age_in_years_non_Z <- as.numeric(dat_no_compt$age_in_years)
dat_no_compt$age_in_years <- scale(as.numeric(dat_no_compt$age_in_years))
dat_no_compt$age_group_coded <- as.factor(dat_no_compt$age_group_coded)
dat_no_compt$RT_info_choice <- scale(as.numeric(dat_no_compt$RT_info_choice))
dat_no_compt$subject_ID <- as.factor(dat_no_compt$subject_ID)
dat_no_compt$gender_coded[dat_no_compt$gender_coded==3] <- 1 # 1 = male, 2 = female, 3 = other
dat_no_compt$gender_coded <- as.factor(dat_no_compt$gender_coded)
contrasts(dat_no_compt$gender_coded) <- contr.helmert(2)
dat_no_compt$info_choice <- as.factor(dat_no_compt$info_choice) # 0 = left, 1 = right
dat_no_compt$percent_comprehension_non_Z <- as.numeric(dat_no_compt$percent_comprehension)
dat_no_compt$percent_comprehension <- scale(as.numeric(dat_no_compt$percent_comprehension))
dat_no_compt$wob_non_Z <- as.numeric(dat_no_compt$wob)
dat_no_compt$wob <- scale(as.numeric(dat_no_compt$wob))
dat_no_compt$chance <- 0.5

## SE1 Compute scaled deltas (chidlren)
# Find min and max and rescale between -1 and 1
dat_no_compt <- rescale_deltas(dat_no_compt)

# Check rescaled and raw deltas have the same sign
check_delta_signs(dat_no_compt)

## SE1 - adults
adu_dat_no_compt <- read.csv("data_supplementary_experiment_1_adults.csv")
adu_dat_no_compt <- subset(adu_dat_no_compt, !(condition %in% c("catch_1", "catch_2")))
adu_dat_no_compt <- subset(adu_dat_no_compt, (catch_trials_score == 100) & (wob > 0.5))

# Turn variables into appropriate data types
adu_dat_no_compt$delta_EV_non_Z <- as.numeric(adu_dat_no_compt$delta_EV)
adu_dat_no_compt$delta_uncertainty_level_non_Z <- as.numeric(adu_dat_no_compt$delta_uncertainty_level)
adu_dat_no_compt$delta_agency_non_Z <- as.numeric(adu_dat_no_compt$delta_agency)
adu_dat_no_compt$age_in_years_non_Z <- as.numeric(adu_dat_no_compt$age_in_years)
adu_dat_no_compt$age_in_years <- scale(as.numeric(adu_dat_no_compt$age_in_years))
adu_dat_no_compt$age_group <- "Adult"
adu_dat_no_compt$age_group_coded <- as.factor(4)
adu_dat_no_compt$RT_info_choice <- scale(as.numeric(adu_dat_no_compt$RT_info_choice))
adu_dat_no_compt$subject_ID <- as.factor(adu_dat_no_compt$subject_ID)
adu_dat_no_compt$gender_coded <- as.factor(adu_dat_no_compt$gender_coded)
contrasts(adu_dat_no_compt$gender_coded) <- contr.helmert(2)
adu_dat_no_compt$info_choice <- as.factor(adu_dat_no_compt$info_choice) # 0 = left, 1 = right
adu_dat_no_compt$percent_comprehension_non_Z <- as.numeric(adu_dat_no_compt$percent_comprehension)
adu_dat_no_compt$percent_comprehension <- scale(as.numeric(adu_dat_no_compt$percent_comprehension))
adu_dat_no_compt$wob_non_Z <- as.numeric(adu_dat_no_compt$wob)
adu_dat_no_compt$wob <- scale(as.numeric(adu_dat_no_compt$wob))
adu_dat_no_compt$chance <- 0.5

## SE1 Compute scaled deltas (chidlren)
# Find min and max and rescale between -1 and 1
adu_dat_no_compt <- rescale_deltas(adu_dat_no_compt)

# Check rescaled and raw deltas have the same sign
check_delta_signs(adu_dat_no_compt)

## SE1 Merged data sets (children + adults) 
cols <- c("subject_ID", "info_choice", 
          "delta_EV", "delta_uncertainty_level", "delta_agency", 
          "EV_L", "EV_R", "uncertainty_level_L", "uncertainty_level_R",
          "agency_probL", "agency_probR",
          "gender_coded", "wob", "percent_comprehension", 
          "wob_non_Z", "percent_comprehension_non_Z", 
          "age_group", "age_group_coded")
dat_no_compt_child_adu <- rbind(dplyr::select(dat_no_compt, all_of(cols)), 
                                dplyr::select(adu_dat_no_compt, all_of(cols)))

# Rescale and re-Z-score variables
dat_no_compt_child_adu$percent_comprehension <- scale(dat_no_compt_child_adu$percent_comprehension_non_Z)
dat_no_compt_child_adu$wob <- scale(dat_no_compt_child_adu$wob_non_Z)
dat_no_compt_child_adu <- rescale_deltas(dat_no_compt_child_adu)

# Fix gender contrasts
contrasts(dat_no_compt_child_adu$gender_coded) <- contr.helmert(2)


# Create new variable and contrasts for children vs adults factor
dat_no_compt_child_adu$group <- ifelse(dat_no_compt_child_adu$age_group_coded == 4, 1, -1)

#### SE1 Demographics ####
summary(subset(dat_no_compt, condition=="affective_1")$age_in_years_non_Z)
sdamr::sample_sd(subset(dat_no_compt, condition=="affective_1")$age_in_years_non_Z)
psych::describeBy(subset(dat_no_compt, condition=="affective_1")$age_in_years_non_Z, 
           subset(dat_no_compt, condition=="affective_1")$age_group)
table(subset(dat_no_compt, condition=="affective_1")$gender)
table(subset(dat_no_compt, condition=="affective_1")$gender, 
      subset(dat_no_compt, condition=="affective_1")$age_group)
psych::describeBy(subset(adu_dat_no_compt, condition=="affective_1")$age_in_years_non_Z, 
           subset(adu_dat_no_compt, condition=="affective_1")$age_group)
table(subset(adu_dat_no_compt, condition=="affective_1")$gender, 
      subset(adu_dat_no_compt, condition=="affective_1")$age_group)

#### SE1 Comprehension and wob scores ####
## t-tests

# Within group
t.test(subset(dat_no_compt, condition == "affective_1")$wob_non_Z, mu=0.5, var.equal=TRUE)
t.test(subset(adu_dat_no_compt, condition == "affective_1")$wob_non_Z, mu=0.5, var.equal=TRUE)

# Between groups
t.test(subset(dat_no_compt, condition == "affective_1")$percent_comprehension_non_Z, 
       subset(adu_dat_no_compt, condition == "affective_1")$percent_comprehension_non_Z, var.equal=TRUE)
t.test(subset(dat_no_compt, condition == "affective_1")$wob_non_Z, 
       subset(adu_dat_no_compt, condition == "affective_1")$wob_non_Z, var.equal=TRUE)

# Standard deviations
sdamr::sample_sd(subset(dat_no_compt, condition == "affective_1")$percent_comprehension_non_Z)
sdamr::sample_sd(subset(dat_no_compt, condition == "affective_1")$wob_non_Z)
sdamr::sample_sd(subset(adu_dat_no_compt, condition == "affective_1")$percent_comprehension_non_Z)
sdamr::sample_sd(subset(adu_dat_no_compt, condition == "affective_1")$wob_non_Z)


## Create long form data sets

# Create separate data sets with variables for children and adults' pc and wob
dat_no_compt_pc_wob <- dplyr::select(subset(dat_no_compt, condition == "affective_1"), 
                                  all_of(c("subject_ID", "age_group_coded", 
                                           "percent_comprehension_non_Z", 
                                           "wob_non_Z")))
dat_no_compt_pc_wob$percent_comprehension_child <- dat_no_compt_pc_wob$"percent_comprehension_non_Z"
dat_no_compt_pc_wob$wob_child <- dat_no_compt_pc_wob$"wob_non_Z"
adu_dat_no_compt_pc_wob <- dplyr::select(subset(adu_dat_no_compt,condition == "affective_1"), 
                                      all_of(c("subject_ID", "age_group_coded", 
                                               "percent_comprehension_non_Z", 
                                               "wob_non_Z")))

# Merge datasets
dat_no_compt_child_adu_pc_wob <- rbind(dplyr::select(dat_no_compt_pc_wob, all_of(
  c("subject_ID", "age_group_coded", "percent_comprehension_non_Z", "wob_non_Z"))), 
  dplyr::select(adu_dat_no_compt_pc_wob, all_of(
    c("subject_ID", "age_group_coded", "percent_comprehension_non_Z", "wob_non_Z"))))

# Create group factor 
dat_no_compt_child_adu_pc_wob$group <- as.factor(ifelse(dat_no_compt_child_adu_pc_wob$age_group_coded == 4, "Adults", "Children"))

# Turn into long form 
dat_no_compt_child_adu_pc_wob <- tidyr::gather(dat_no_compt_child_adu_pc_wob, key = "variable", value = "score", percent_comprehension_non_Z:wob_non_Z, factor_key=TRUE)

# Create a new variable to distinguish between children's and adult's scores
dat_no_compt_child_adu_pc_wob$variable2 <- "None"
dat_no_compt_child_adu_pc_wob$variable2[dat_no_compt_child_adu_pc_wob$group == "Children" 
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_child"
dat_no_compt_child_adu_pc_wob$variable2[dat_no_compt_child_adu_pc_wob$group == "Children" 
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_child"
dat_no_compt_child_adu_pc_wob$variable2[dat_no_compt_child_adu_pc_wob$group == "Adults" 
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_adu"
dat_no_compt_child_adu_pc_wob$variable2[dat_no_compt_child_adu_pc_wob$group == "Adults" 
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_adu"
dat_no_compt_child_adu_pc_wob$variable2 <- factor(dat_no_compt_child_adu_pc_wob$variable2, 
                                               levels=c("percent_comprehension_child",
                                                        "percent_comprehension_adu",
                                                        "wob_child",
                                                        "wob_adu"))
# Create a new variable for colors for variable2
dat_no_compt_child_adu_pc_wob$color <- "None"
dat_no_compt_child_adu_pc_wob$color <- as.factor(ifelse(
  dat_no_compt_child_adu_pc_wob$age_group_coded == 4, "#444444", "#AAAAAA"))



# Create a new variable to distinguish between each group of children's and adult's scores
dat_no_compt_child_adu_pc_wob$variable3 <- "None"
dat_no_compt_child_adu_pc_wob$variable3[dat_no_compt_child_adu_pc_wob$age_group_coded == 0
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_4_5"
dat_no_compt_child_adu_pc_wob$variable3[dat_no_compt_child_adu_pc_wob$age_group_coded == 0
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_4_5"
dat_no_compt_child_adu_pc_wob$variable3[dat_no_compt_child_adu_pc_wob$age_group_coded == 1
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_6_7"
dat_no_compt_child_adu_pc_wob$variable3[dat_no_compt_child_adu_pc_wob$age_group_coded == 1
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_6_7"
dat_no_compt_child_adu_pc_wob$variable3[dat_no_compt_child_adu_pc_wob$age_group_coded == 2
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_8_9"
dat_no_compt_child_adu_pc_wob$variable3[dat_no_compt_child_adu_pc_wob$age_group_coded == 2
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_8_9"
dat_no_compt_child_adu_pc_wob$variable3[dat_no_compt_child_adu_pc_wob$age_group_coded == 3
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_10_12"
dat_no_compt_child_adu_pc_wob$variable3[dat_no_compt_child_adu_pc_wob$age_group_coded == 3
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_10_12"
dat_no_compt_child_adu_pc_wob$variable3[dat_no_compt_child_adu_pc_wob$age_group_coded == 4 
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "percent_comprehension_non_Z"] <- "percent_comprehension_adu"
dat_no_compt_child_adu_pc_wob$variable3[dat_no_compt_child_adu_pc_wob$age_group_coded == 4
                                     & dat_no_compt_child_adu_pc_wob$variable ==
                                       "wob_non_Z"] <- "wob_adu"
dat_no_compt_child_adu_pc_wob$variable3 <- factor(dat_no_compt_child_adu_pc_wob$variable3, 
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

table(dat_no_compt_child_adu_pc_wob$variable3)

## Plots
# Children vs adults
no_compt_pc <- ggplot(subset(dat_no_compt_child_adu_pc_wob, variable2 %in% 
                            c("percent_comprehension_child", "percent_comprehension_adu")), 
                   aes(x=variable2, y=score, fill=variable2)) +
  geom_jitter(height=0, color = subset(dat_no_compt_child_adu_pc_wob, 
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
no_compt_pc


no_compt_wob <- ggplot(subset(dat_no_compt_child_adu_pc_wob, variable2 %in% 
                             c("wob_child", "wob_adu")), 
                    aes(x=variable2, y=score, fill=variable2)) +
  geom_jitter(height=0, color = subset(dat_no_compt_child_adu_pc_wob, 
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
no_compt_wob

grid.arrange(no_compt_pc, no_compt_wob, ncol=2)

# By group
no_compt_pc_by_group <- ggplot(subset(dat_no_compt_child_adu_pc_wob, variable2 %in% 
                                     c("percent_comprehension_child", "percent_comprehension_adu")), 
                            aes(x=variable3, y=score, fill=variable2)) +
  geom_jitter(height=0, color = subset(dat_no_compt_child_adu_pc_wob, 
                                       variable2 %in% c("percent_comprehension_child", "percent_comprehension_adu"))$color,
              size=1.5, alpha=0.9)+   
  geom_violin(alpha=0.3) +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0, 1)) + 
  scale_fill_manual(values=c("#AAAAAA", "#444444")) +
  #stat_summary(fun=mean, geom="point", shape=20, size=4, color="black", fill="black") +
  theme_classic() + labs(x="comprehension", y="Score") + theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm"))
no_compt_pc_by_group

no_compt_wob_by_group <- ggplot(subset(dat_no_compt_child_adu_pc_wob, variable2 %in% 
                                      c("wob_child", "wob_adu")), 
                             aes(x=variable3, y=score, fill=variable2)) +
  geom_jitter(height=0, color = subset(dat_no_compt_child_adu_pc_wob, 
                                       variable2 %in% c("wob_child", "wob_adu"))$color,
              size=1.5, alpha=0.9)+   
  geom_violin(alpha=0.3) +
  geom_hline(yintercept=0.5, color="#AAAAAA") +
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0, 1)) + 
  scale_fill_manual(values=c("#AAAAAA", "#444444")) +
  #stat_summary(fun=mean, geom="point", shape=20, size=4, color="black", fill="black") +
  theme_classic() + labs(x="comprehension", y="Score") + theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm"))
no_compt_wob_by_group

grid.arrange(no_compt_pc_by_group, no_compt_wob_by_group, ncol=2)


#### SE1 Models by age group ####
# Children overall 
# Chance only
child_no_compt_mod_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID),
                                data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, with interactions
child_no_compt_mod_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                percent_comprehension + wob + gender_coded + age_in_years +
                                age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                (delta_EV + delta_uncertainty_level + delta_agency
                                 | subject_ID), 
                              data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, without intercept
child_no_compt_mod_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                     percent_comprehension + wob + gender_coded + age_in_years - 1 +
                                                     age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                                     (delta_EV + delta_uncertainty_level + delta_agency
                                                      | subject_ID), 
                                                   data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full,  without interaction between delta_agency and age
child_no_compt_mod_full_drop_Af_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                            percent_comprehension + wob + gender_coded + age_in_years +
                                            age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                            (delta_EV + delta_uncertainty_level + delta_agency
                                             | subject_ID), 
                                          data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without interaction between delta_uncertainty_level and age
child_no_compt_mod_full_drop_C_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                           percent_comprehension + wob + gender_coded + age_in_years +
                                           age_in_years:delta_EV + age_in_years:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency
                                            | subject_ID), 
                                         data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without interaction between delta_agency and age
child_no_compt_mod_full_drop_Ac_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                            percent_comprehension + wob + gender_coded + age_in_years +
                                            age_in_years:delta_EV + age_in_years:delta_uncertainty_level +
                                            (delta_EV + delta_uncertainty_level + delta_agency
                                             | subject_ID), 
                                          data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without EV as a fixed effect
child_no_compt_mod_full_drop_Af_fixed <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                                              percent_comprehension + wob + gender_coded + age_in_years +
                                              age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without uncertainty as a fixed effect
child_no_compt_mod_full_drop_C_fixed <- glmer(info_choice ~ delta_EV + delta_agency +
                                             percent_comprehension + wob + gender_coded + age_in_years +
                                             age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                             (delta_EV + delta_uncertainty_level + delta_agency
                                              | subject_ID), 
                                           data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without agency as a fixed effect
child_no_compt_mod_full_drop_Ac_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                              percent_comprehension + wob + gender_coded + age_in_years +
                                              age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without comprehension as a fixed effect
child_no_compt_mod_full_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                              wob + gender_coded + age_in_years +
                                              age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                              (delta_EV + delta_uncertainty_level + delta_agency
                                               | subject_ID), 
                                            data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without wob as a fixed effect
child_no_compt_mod_full_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + gender_coded + age_in_years +
                                               age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                               (delta_EV + delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without gender as a fixed effect
child_no_compt_mod_full_drop_gender_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                  percent_comprehension + wob  + age_in_years +
                                                  age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                                  (delta_EV + delta_uncertainty_level + delta_agency
                                                   | subject_ID), 
                                                data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, without age as a fixed effect
child_no_compt_mod_full_drop_age_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                               percent_comprehension + wob + gender_coded +
                                               age_in_years:delta_EV + age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                               (delta_EV + delta_uncertainty_level + delta_agency
                                                | subject_ID), 
                                             data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, without EV as a fixed or interaction effect
child_no_compt_mod_full_drop_Af_fixed_and_int <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                                                      percent_comprehension + wob + gender_coded + age_in_years +
                                                      age_in_years:delta_uncertainty_level + age_in_years:delta_agency +
                                                      (delta_EV + delta_uncertainty_level + delta_agency
                                                       | subject_ID), 
                                                    data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without EV as a fixed or interaction effect
child_no_compt_mod_full_drop_C_fixed_and_int <- glmer(info_choice ~ delta_EV + delta_agency +
                                                     percent_comprehension + wob + gender_coded + age_in_years +
                                                     age_in_years:delta_EV + age_in_years:delta_agency +
                                                     (delta_EV + delta_uncertainty_level + delta_agency
                                                      | subject_ID), 
                                                   data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without agency as a fixed or interaction effect
child_no_compt_mod_full_drop_Ac_fixed_and_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                                      percent_comprehension + wob + gender_coded + age_in_years +
                                                      age_in_years:delta_EV + age_in_years:delta_uncertainty_level +
                                                      (delta_EV + delta_uncertainty_level + delta_agency
                                                       | subject_ID), 
                                                    data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Affect only
child_no_compt_mod_Af <- glmer(info_choice ~ delta_EV +
                              percent_comprehension + wob + gender_coded + age_in_years +
                              delta_EV:age_in_years +
                              (delta_EV 
                               | subject_ID), 
                            data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition only 
child_no_compt_mod_C <- glmer(info_choice ~  delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + age_in_years +
                             delta_uncertainty_level:age_in_years +
                             (delta_uncertainty_level 
                              | subject_ID),
                           data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Action only
child_no_compt_mod_Ac <- glmer(info_choice ~ delta_agency +
                              percent_comprehension + wob + gender_coded + age_in_years +
                              delta_agency:age_in_years +
                              (delta_agency 
                               | subject_ID), 
                            data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Cognition 
child_no_compt_mod_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                               percent_comprehension + wob + gender_coded + age_in_years +
                               delta_EV:age_in_years + delta_uncertainty_level:age_in_years +
                               (delta_EV + delta_uncertainty_level 
                                | subject_ID),
                             data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Action 
child_no_compt_mod_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                                percent_comprehension + wob + gender_coded + age_in_years +
                                delta_EV:age_in_years + delta_agency:age_in_years +
                                (delta_EV + delta_agency
                                 | subject_ID),
                              data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition +  Action
child_no_compt_mod_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + age_in_years +
                               delta_uncertainty_level:age_in_years + delta_agency:age_in_years +
                               (delta_uncertainty_level + delta_agency 
                                | subject_ID),
                             data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Adults
# Chance only
adu_no_compt_mod_chance <- glmer(info_choice ~ 0 + chance + (1|subject_ID), 
                              data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full 
adu_no_compt_mod_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_uncertainty_level + delta_agency 
                               | subject_ID), 
                            data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, without intercept
adu_no_compt_mod_full_drop_intercept_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                   percent_comprehension + wob + gender_coded -1 +
                                                   (delta_EV + delta_uncertainty_level + delta_agency 
                                                    | subject_ID), 
                                                 data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, without EV as a fixed effect
adu_no_compt_mod_full_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                            percent_comprehension + wob + gender_coded + 
                                            (delta_EV + delta_uncertainty_level + delta_agency 
                                             | subject_ID), 
                                          data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without uncertainty as a fixed effect
adu_no_compt_mod_full_drop_C_fixed <- glmer(info_choice ~  delta_EV + delta_agency +
                                           percent_comprehension + wob + gender_coded + 
                                           (delta_EV + delta_uncertainty_level + delta_agency 
                                            | subject_ID), 
                                         data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without agency as a fixed effect
adu_no_compt_mod_full_drop_Ac_fixed <- glmer(info_choice ~  delta_EV + delta_uncertainty_level +
                                            percent_comprehension + wob + gender_coded + 
                                            (delta_EV + delta_uncertainty_level + delta_agency 
                                             | subject_ID), 
                                          data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Full, without comprehension as a fixed effect
adu_no_compt_mod_full_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                            wob + gender_coded +
                                            (delta_EV + delta_uncertainty_level + delta_agency
                                             | subject_ID), 
                                          data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without wob as a fixed effect
adu_no_compt_mod_full_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                             percent_comprehension + gender_coded +
                                             (delta_EV + delta_uncertainty_level + delta_agency
                                              | subject_ID), 
                                           data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Full, without gender as a fixed effect
adu_no_compt_mod_full_drop_gender_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                                percent_comprehension + wob  +
                                                (delta_EV + delta_uncertainty_level + delta_agency
                                                 | subject_ID), 
                                              data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)


# Affect only
adu_no_compt_mod_Af <- glmer(info_choice ~ delta_EV +
                            percent_comprehension + wob + gender_coded + 
                            (delta_EV 
                             | subject_ID), 
                          data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition only
adu_no_compt_mod_C <- glmer(info_choice ~  delta_uncertainty_level +
                           percent_comprehension + wob + gender_coded + 
                           (delta_uncertainty_level 
                            | subject_ID), 
                         data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Action only 
adu_no_compt_mod_Ac <- glmer(info_choice ~ delta_agency +
                            percent_comprehension + wob + gender_coded +
                            (delta_agency 
                             | subject_ID), 
                          data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Cognition 
adu_no_compt_mod_AfC <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                             percent_comprehension + wob + gender_coded + 
                             (delta_EV + delta_uncertainty_level 
                              | subject_ID), 
                           data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Affect + Action 
adu_no_compt_mod_AfAc <- glmer(info_choice ~ delta_EV + delta_agency +
                              percent_comprehension + wob + gender_coded + 
                              (delta_EV + delta_agency
                               | subject_ID), 
                            data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Cognition + Action 
adu_no_compt_mod_CAc <- glmer(info_choice ~ delta_uncertainty_level + delta_agency +
                             percent_comprehension + wob + gender_coded + 
                             (delta_uncertainty_level + delta_agency 
                              | subject_ID), 
                           data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

#### SE1 Significance tests within groups ####
# (within single groups)
# access results of test a with a$Chisq, a$Df, a$Pr[2]

## Children overall
# delta EV: overall, only fixed, fixed + interaction
anova(child_no_compt_mod_full, child_no_compt_mod_CAc, test="Chisq")
child_test_fixed_EV <- anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_Af_fixed, test="Chisq")
anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_Af_fixed_and_int, test="Chisq")
# delta uncertainty: overall, only fixed, fixed + interaction
anova(child_no_compt_mod_full, child_no_compt_mod_AfAc, test="Chisq")
child_test_fixed_uncertainty <- anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_C_fixed, test="Chisq")
anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_C_fixed_and_int, test="Chisq")
# delta agency: overall, only fixed, fixed + interaction
anova(child_no_compt_mod_full, child_no_compt_mod_AfC, test="Chisq")
child_test_fixed_agency <-anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_Ac_fixed, test="Chisq")
anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_Ac_fixed_and_int, test="Chisq")
# delta EV x age
child_test_fixed_EV_int <-anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_Af_int, test="Chisq")
# delta uncertainty x age
child_test_fixed_uncertainty_int <- anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_C_int, test="Chisq")
# delta agency x age
child_test_fixed_agency_int <- anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_Ac_int, test="Chisq")
# covariates
child_test_fixed_pc <- anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_pc_fixed, test="Chisq")
child_test_fixed_wob <- anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_wob_fixed, test="Chisq")
child_test_fixed_gender <- anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_gender_fixed, test="Chisq")
child_test_fixed_age <- anova(child_no_compt_mod_full, child_no_compt_mod_full_drop_age_fixed, test="Chisq")

## Adults overall
# delta EV: overall, only fixed
anova(adu_no_compt_mod_full, adu_no_compt_mod_CAc, test="Chisq")
adu_test_fixed_EV <- anova(adu_no_compt_mod_full, adu_no_compt_mod_full_drop_Af_fixed, test="Chisq")
# delta uncertainty: overall, only fixed
anova(adu_no_compt_mod_full, adu_no_compt_mod_AfAc, test="Chisq")
adu_test_fixed_uncertainty <- anova(adu_no_compt_mod_full, adu_no_compt_mod_full_drop_C_fixed, test="Chisq")
# delta agency: overall, only fixed
anova(adu_no_compt_mod_full, adu_no_compt_mod_AfC, test="Chisq")
adu_test_fixed_agency <- anova(adu_no_compt_mod_full, adu_no_compt_mod_full_drop_Ac_fixed, test="Chisq")
# covariates
adu_test_fixed_pc <- anova(adu_no_compt_mod_full, adu_no_compt_mod_full_drop_pc_fixed, test="Chisq")
adu_test_fixed_wob <- anova(adu_no_compt_mod_full, adu_no_compt_mod_full_drop_wob_fixed, test="Chisq")
adu_test_fixed_gender <- anova(adu_no_compt_mod_full, adu_no_compt_mod_full_drop_gender_fixed, test="Chisq")


# Create csv file with tests of fixed effects
variable <- c("children_EV",  "children_uncertainty", "children_agency", 
              "children_comprehension", "children_wob", "children_gender", 
              "children_age", "children_age_by_EV", 
              "children_age_by_uncertainty", "children_age_by_agency",
              "adults_EV",  "adults_uncertainty", "adults_agency", 
              "adults_comprehension", "adults_wob", "adults_gender")

# store coefficients
child_no_compt_mod_full_coeffs <- summary(child_no_compt_mod_full)$coefficients[, 1]
adu_no_compt_mod_full_coeffs <- summary(adu_no_compt_mod_full)$coefficients[, 1]

# store sems
child_no_compt_mod_full_sems <- summary(child_no_compt_mod_full)$coefficients[, 2]
adu_no_compt_mod_full_sems <- summary(adu_no_compt_mod_full)$coefficients[, 2]

estimate <- c(child_no_compt_mod_full_coeffs["delta_EV"], 
              child_no_compt_mod_full_coeffs["delta_uncertainty_level"],
              child_no_compt_mod_full_coeffs["delta_agency"],
              child_no_compt_mod_full_coeffs["percent_comprehension"],
              child_no_compt_mod_full_coeffs["wob"],
              child_no_compt_mod_full_coeffs["gender_coded1"],
              child_no_compt_mod_full_coeffs["age_in_years"],
              child_no_compt_mod_full_coeffs["delta_EV:age_in_years"],
              child_no_compt_mod_full_coeffs["delta_uncertainty_level:age_in_years"],
              child_no_compt_mod_full_coeffs["delta_agency:age_in_years"],
              adu_no_compt_mod_full_coeffs["delta_EV"], 
              adu_no_compt_mod_full_coeffs["delta_uncertainty_level"],
              adu_no_compt_mod_full_coeffs["delta_agency"],
              adu_no_compt_mod_full_coeffs["percent_comprehension"],
              adu_no_compt_mod_full_coeffs["wob"],
              adu_no_compt_mod_full_coeffs["gender_coded1"])

sem <- c(child_no_compt_mod_full_sems["delta_EV"], 
         child_no_compt_mod_full_sems["delta_uncertainty_level"],
         child_no_compt_mod_full_sems["delta_agency"],
         child_no_compt_mod_full_sems["percent_comprehension"],
         child_no_compt_mod_full_sems["wob"],
         child_no_compt_mod_full_sems["gender_coded1"],
         child_no_compt_mod_full_sems["age_in_years"],
         child_no_compt_mod_full_sems["delta_EV:age_in_years"],
         child_no_compt_mod_full_sems["delta_uncertainty_level:age_in_years"],
         child_no_compt_mod_full_sems["delta_agency:age_in_years"],
         adu_no_compt_mod_full_sems["delta_EV"], 
         adu_no_compt_mod_full_sems["delta_uncertainty_level"],
         adu_no_compt_mod_full_sems["delta_agency"],
         adu_no_compt_mod_full_sems["percent_comprehension"],
         adu_no_compt_mod_full_sems["wob"],
         adu_no_compt_mod_full_sems["gender_coded1"])


chisq <- c(child_test_fixed_EV$Chisq[2], child_test_fixed_uncertainty$Chisq[2],
           child_test_fixed_agency$Chisq[2], child_test_fixed_pc$Chisq[2],
           child_test_fixed_wob$Chisq[2], child_test_fixed_gender$Chisq[2], 
           child_test_fixed_age$Chisq[2], child_test_fixed_EV_int$Chisq[2], 
           child_test_fixed_uncertainty_int$Chisq[2], child_test_fixed_agency_int$Chisq[2],
           adu_test_fixed_EV$Chisq[2], adu_test_fixed_uncertainty$Chisq[2],
           adu_test_fixed_agency$Chisq[2], adu_test_fixed_pc$Chisq[2],
           adu_test_fixed_wob$Chisq[2], adu_test_fixed_gender$Chisq[2])

df <- c(child_test_fixed_EV$Df[2], child_test_fixed_uncertainty$Df[2],
        child_test_fixed_agency$Df[2], child_test_fixed_pc$Df[2],
        child_test_fixed_wob$Df[2], child_test_fixed_gender$Df[2], 
        child_test_fixed_age$Df[2], child_test_fixed_EV_int$Df[2], 
        child_test_fixed_uncertainty_int$Df[2], child_test_fixed_agency_int$Df[2],
        adu_test_fixed_EV$Df[2], adu_test_fixed_uncertainty$Df[2],
        adu_test_fixed_agency$Df[2], adu_test_fixed_pc$Df[2],
        adu_test_fixed_wob$Df[2], adu_test_fixed_gender$Df[2])

p_value <- c(child_test_fixed_EV$Pr[2], child_test_fixed_uncertainty$Pr[2],
             child_test_fixed_agency$Pr[2], child_test_fixed_pc$Pr[2],
             child_test_fixed_wob$Pr[2], child_test_fixed_gender$Pr[2], 
             child_test_fixed_age$Pr[2], child_test_fixed_EV_int$Pr[2], 
             child_test_fixed_uncertainty_int$Pr[2], child_test_fixed_agency_int$Pr[2],
             adu_test_fixed_EV$Pr[2], adu_test_fixed_uncertainty$Pr[2],
             adu_test_fixed_agency$Pr[2], adu_test_fixed_pc$Pr[2],
             adu_test_fixed_wob$Pr[2], adu_test_fixed_gender$Pr[2])


no_compt_mods_fixed_effects <- data.frame(variable, estimate, sem, chisq, df, p_value)

# Save fixed effects as a file
write.csv(no_compt_mods_fixed_effects, sprintf("%s/no_compt_mods_fixed_effects.csv", wd), row.names = FALSE)

#### SE1 AIC ####
# Children vs adults
age_group <- factor(c("Children", "Adults"), levels=c("Children", "Adults"))
chance_AIC <- c(AIC(child_no_compt_mod_chance), AIC(adu_no_compt_mod_chance))
Af_AIC <- c(AIC(child_no_compt_mod_Af), AIC(adu_no_compt_mod_Af))
C_AIC <- c(AIC(child_no_compt_mod_C), AIC(adu_no_compt_mod_C))
Ac_AIC <- c(AIC(child_no_compt_mod_Ac), AIC(adu_no_compt_mod_Ac))
AfC_AIC <- c(AIC(child_no_compt_mod_AfC), AIC(adu_no_compt_mod_AfC))
AfAc_AIC <- c(AIC(child_no_compt_mod_AfAc), AIC(adu_no_compt_mod_AfAc))
CAc_AIC <- c(AIC(child_no_compt_mod_CAc), AIC(adu_no_compt_mod_CAc))
full_AIC <- c(AIC(child_no_compt_mod_full), AIC(adu_no_compt_mod_full))

AIC_scores_no_compt_child_vs_adu <- data.frame(age_group, chance_AIC, 
                                            Af_AIC, C_AIC, Ac_AIC, AfC_AIC, AfAc_AIC, CAc_AIC, full_AIC)
# With covariates and random effects
write.csv(AIC_scores_no_compt_child_vs_adu, sprintf("%s/AIC_scores_no_compt_child_vs_adu.csv", wd), row.names = FALSE)

AIC_scores_no_compt_child_vs_adu_long <- gather(AIC_scores_no_compt_child_vs_adu, model, AIC_score, chance_AIC:full_AIC, factor_key=TRUE)
AIC_scores_no_compt_child_vs_adu_long$model <- factor(AIC_scores_no_compt_child_vs_adu_long$model, labels=c("Chance", 
                                                                                                      "EV", "Uncertainty", "Agency",
                                                                                                      "EV + Uncertainty", "EV + Agency", "Uncertainty + Agency", "EV + Uncertainty + Agency"))
AIC_scores_no_compt_child_vs_adu_plot <- ggplot(data=AIC_scores_no_compt_child_vs_adu_long, aes(x=model, y=AIC_score)) +
  geom_bar(stat="identity") + 
  facet_grid(~age_group) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=60, hjust=1, size=14)) +
  theme(axis.text.y = element_text(size=14))
AIC_scores_no_compt_child_vs_adu_plot


#### SE1 Betas plots #### 
# Create betas plots
# Children vs adults
no_compt_hedonic_beta <- c(summary(child_no_compt_mod_full)$coefficients[2], summary(adu_no_compt_mod_full)$coefficients[2])
no_compt_cognitive_beta <- c(summary(child_no_compt_mod_full)$coefficients[3], summary(adu_no_compt_mod_full)$coefficients[3])
no_compt_instrumental_beta <- c(summary(child_no_compt_mod_full)$coefficients[4], summary(adu_no_compt_mod_full)$coefficients[4])
no_compt_hedonic_sem_lower <- c(summary(child_no_compt_mod_full)$coefficients[2] - summary(child_no_compt_mod_full)$coefficients[2,2],
                             summary(adu_no_compt_mod_full)$coefficients[2] - summary(adu_no_compt_mod_full)$coefficients[2,2])
no_compt_cognitive_sem_lower <- c(summary(child_no_compt_mod_full)$coefficients[3] - summary(child_no_compt_mod_full)$coefficients[3,2],
                               summary(adu_no_compt_mod_full)$coefficients[3] - summary(adu_no_compt_mod_full)$coefficients[3,2])
no_compt_instrumental_sem_lower <- c(summary(child_no_compt_mod_full)$coefficients[4] - summary(child_no_compt_mod_full)$coefficients[4,2],
                                  summary(adu_no_compt_mod_full)$coefficients[4] - summary(adu_no_compt_mod_full)$coefficients[4,2])
no_compt_hedonic_sem_upper <- c(summary(child_no_compt_mod_full)$coefficients[2] + summary(child_no_compt_mod_full)$coefficients[2,2],
                             summary(adu_no_compt_mod_full)$coefficients[2] + summary(adu_no_compt_mod_full)$coefficients[2,2])
no_compt_cognitive_sem_upper <- c(summary(child_no_compt_mod_full)$coefficients[3] + summary(child_no_compt_mod_full)$coefficients[3,2],
                               summary(adu_no_compt_mod_full)$coefficients[3] + summary(adu_no_compt_mod_full)$coefficients[3,2])
no_compt_instrumental_sem_upper <- c(summary(child_no_compt_mod_full)$coefficients[4] + summary(child_no_compt_mod_full)$coefficients[4,2],
                                  summary(adu_no_compt_mod_full)$coefficients[4] + summary(adu_no_compt_mod_full)$coefficients[4,2])
hedonic_color <- "#E41A1C"
cognitive_color <- "#377EB8"
instrumental_color <- "#4DAF4A"

no_compt_betas_dat <- data.frame(no_compt_hedonic_beta, no_compt_cognitive_beta, no_compt_instrumental_beta,
                              no_compt_hedonic_sem_lower, no_compt_cognitive_sem_lower, no_compt_instrumental_sem_lower,
                              no_compt_hedonic_sem_upper, no_compt_cognitive_sem_upper, no_compt_instrumental_sem_upper,
                              hedonic_color, cognitive_color, instrumental_color)

no_compt_motive_betas <- ggplot()+
  scale_y_continuous(breaks=seq(0,3.5,0.1), limits=c(0, 3.5), expand=c(0, 0)) + 
  # Children
  geom_errorbar(aes(0.3, ymin = no_compt_hedonic_sem_lower[1], ymax = no_compt_hedonic_sem_upper[1]), size = 0.8, width = .05, colour = hedonic_color) +
  geom_point(aes(0.3, no_compt_hedonic_beta[1]), size = 4, color = hedonic_color)+
  geom_errorbar(aes(0.6, ymin = no_compt_cognitive_sem_lower[1], ymax = no_compt_cognitive_sem_upper[1]), size = 0.8, width = .05, colour = cognitive_color) +
  geom_point(aes(0.6, no_compt_cognitive_beta[1]), size = 4, color = cognitive_color)+
  geom_errorbar(aes(0.9, ymin = no_compt_instrumental_sem_lower[1], ymax = no_compt_instrumental_sem_upper[1]), size = 0.8, width = .05, colour = instrumental_color) +
  geom_point(aes(0.9, no_compt_instrumental_beta[1]), size = 4, color = instrumental_color)+
  # Adults
  geom_errorbar(aes(0.4, ymin = no_compt_hedonic_sem_lower[2], ymax = no_compt_hedonic_sem_upper[2]), size = 0.8, width = .05, colour = hedonic_color, linetype="dashed") +
  geom_point(aes(0.4, no_compt_hedonic_beta[2]), size = 4, color = hedonic_color)+
  geom_errorbar(aes(0.7, ymin = no_compt_cognitive_sem_lower[2], ymax = no_compt_cognitive_sem_upper[2]), size = 0.8, width = .05, colour = cognitive_color, linetype="dashed") +
  geom_point(aes(0.7, no_compt_cognitive_beta[2]), size = 4, color = cognitive_color)+
  geom_errorbar(aes(1, ymin = no_compt_instrumental_sem_lower[2], ymax = no_compt_instrumental_sem_upper[2]), size = 0.8, width = .05, colour = instrumental_color, linetype="dashed") +
  geom_point(aes(1, no_compt_instrumental_beta[2]), size = 4, color = instrumental_color)+
  geom_vline(xintercept=c(0.5, 0.8), linetype="dotted") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size=20)) +
  labs(x="Utilities", y="Standardized beta predicting information-seeking")
no_compt_motive_betas

#### SE1 Children vs adults tests ####
## Models
no_compt_mod_child_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                               percent_comprehension + wob + gender_coded + group + 
                               group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                               (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                             data = dat_no_compt_child_adu, family = binomial, 
                             control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between EV and group
no_compt_mod_child_adu_drop_Af_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                           percent_comprehension + wob + gender_coded + group + 
                                           group:delta_uncertainty_level + group:delta_agency +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_no_compt_child_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between uncertainty and group
no_compt_mod_child_adu_drop_C_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                          percent_comprehension + wob + gender_coded + group + 
                                          group:delta_EV + group:delta_agency +
                                          (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                        data = dat_no_compt_child_adu, family = binomial, 
                                        control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop interaction between agency and group
no_compt_mod_child_adu_drop_Ac_int <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                           percent_comprehension + wob + gender_coded + group + 
                                           group:delta_EV + group:delta_uncertainty_level +
                                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                         data = dat_no_compt_child_adu, family = binomial, 
                                         control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop EV fixed
no_compt_mod_child_adu_drop_Af_fixed <- glmer(info_choice ~  delta_uncertainty_level + delta_agency +
                                  percent_comprehension + wob + gender_coded + group + 
                                  group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                  (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                data = dat_no_compt_child_adu, family = binomial, 
                                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop uncertainty fixed
no_compt_mod_child_adu_drop_C_fixed <- glmer(info_choice ~ delta_EV + delta_agency +
                                  percent_comprehension + wob + gender_coded + group + 
                                  group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                  (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                data = dat_no_compt_child_adu, family = binomial, 
                                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop agency fixed
no_compt_mod_child_adu_drop_Ac_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level +
                                  percent_comprehension + wob + gender_coded + group + 
                                  group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                  (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                data = dat_no_compt_child_adu, family = binomial, 
                                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop pc fixed
no_compt_mod_child_adu_drop_pc_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                   wob + gender_coded + group + 
                                  group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                  (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                data = dat_no_compt_child_adu, family = binomial, 
                                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Drop wob fixed
no_compt_mod_child_adu_drop_wob_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                  percent_comprehension + gender_coded + group + 
                                  group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                  (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                data = dat_no_compt_child_adu, family = binomial, 
                                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop gender fixed
no_compt_mod_child_adu_drop_gender_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                  percent_comprehension + wob + group + 
                                  group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                  (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                data = dat_no_compt_child_adu, family = binomial, 
                                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
# Drop group fixed 
no_compt_mod_child_adu_drop_group_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                  percent_comprehension + wob + gender_coded + 
                                  group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                  (delta_EV + delta_uncertainty_level + delta_agency | subject_ID), 
                                data = dat_no_compt_child_adu, family = binomial, 
                                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)


# Model results
sjPlot::tab_model(no_compt_mod_child_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)

## Test insteractions
# Children versus adults
child_adu_test_EV_int <-anova(no_compt_mod_child_adu, no_compt_mod_child_adu_drop_Af_int, test="Chisq")
child_adu_test_uncertainty_int <-anova(no_compt_mod_child_adu, no_compt_mod_child_adu_drop_C_int, test="Chisq")
child_adu_test_agency_int <-anova(no_compt_mod_child_adu, no_compt_mod_child_adu_drop_Ac_int, test="Chisq")
child_adu_test_fixed_EV <- anova(no_compt_mod_child_adu, no_compt_mod_child_adu_drop_Af_fixed, test="Chisq")
child_adu_test_fixed_uncertainty <- anova(no_compt_mod_child_adu, no_compt_mod_child_adu_drop_C_fixed, test="Chisq")
child_adu_test_fixed_agency <- anova(no_compt_mod_child_adu, no_compt_mod_child_adu_drop_Ac_fixed, test="Chisq")
child_adu_test_fixed_pc <- anova(no_compt_mod_child_adu, no_compt_mod_child_adu_drop_pc_fixed, test="Chisq")
child_adu_test_fixed_wob <- anova(no_compt_mod_child_adu, no_compt_mod_child_adu_drop_wob_fixed, test="Chisq")
child_adu_test_fixed_gender <- anova(no_compt_mod_child_adu, no_compt_mod_child_adu_drop_gender_fixed, test="Chisq")
child_adu_test_fixed_group <- anova(no_compt_mod_child_adu, no_compt_mod_child_adu_drop_group_fixed, test="Chisq")


# Create csv file with tests of fixed effects
variable <- c("child_adu_EV",  "child_adu_uncertainty", "child_adu_agency", 
              "child_adu_comprehension", "child_adu_wob", "child_adu_gender", 
              "child_adu_group", "child_adu_group_by_EV", 
              "child_adu_group_by_uncertainty", "child_adu_group_by_agency")

# store coefficients
child_adu_no_compt_mod_full_coeffs <- summary(no_compt_mod_child_adu)$coefficients[, 1]

# store sems
child_adu_no_compt_mod_full_sems <- summary(no_compt_mod_child_adu)$coefficients[, 2]

estimate <- c(child_adu_no_compt_mod_full_coeffs["delta_EV"], 
              child_adu_no_compt_mod_full_coeffs["delta_uncertainty_level"],
              child_adu_no_compt_mod_full_coeffs["delta_agency"],
              child_adu_no_compt_mod_full_coeffs["percent_comprehension"],
              child_adu_no_compt_mod_full_coeffs["wob"],
              child_adu_no_compt_mod_full_coeffs["gender_coded1"],
              child_adu_no_compt_mod_full_coeffs["group"],
              child_adu_no_compt_mod_full_coeffs["delta_EV:group"],
              child_adu_no_compt_mod_full_coeffs["delta_uncertainty_level:group"],
              child_adu_no_compt_mod_full_coeffs["delta_agency:group"])

sem <- c(child_adu_no_compt_mod_full_sems["delta_EV"], 
              child_adu_no_compt_mod_full_sems["delta_uncertainty_level"],
              child_adu_no_compt_mod_full_sems["delta_agency"],
              child_adu_no_compt_mod_full_sems["percent_comprehension"],
              child_adu_no_compt_mod_full_sems["wob"],
              child_adu_no_compt_mod_full_sems["gender_coded1"],
              child_adu_no_compt_mod_full_sems["group"],
              child_adu_no_compt_mod_full_sems["delta_EV:group"],
              child_adu_no_compt_mod_full_sems["delta_uncertainty_level:group"],
              child_adu_no_compt_mod_full_sems["delta_agency:group"])


chisq <- c(child_adu_test_fixed_EV$Chisq[2], child_adu_test_fixed_uncertainty$Chisq[2],
           child_adu_test_fixed_agency$Chisq[2], child_adu_test_fixed_pc$Chisq[2],
           child_adu_test_fixed_wob$Chisq[2], child_adu_test_fixed_gender$Chisq[2], 
           child_adu_test_fixed_group$Chisq[2], child_adu_test_EV_int$Chisq[2], 
           child_adu_test_uncertainty_int$Chisq[2], child_adu_test_agency_int$Chisq[2])

df <- c(child_adu_test_fixed_EV$Df[2], child_adu_test_fixed_uncertainty$Df[2],
           child_adu_test_fixed_agency$Df[2], child_adu_test_fixed_pc$Df[2],
           child_adu_test_fixed_wob$Df[2], child_adu_test_fixed_gender$Df[2], 
           child_adu_test_fixed_group$Df[2], child_adu_test_EV_int$Df[2], 
           child_adu_test_uncertainty_int$Df[2], child_adu_test_agency_int$Df[2])

p_value <- c(child_adu_test_fixed_EV$Pr[2], child_adu_test_fixed_uncertainty$Pr[2],
           child_adu_test_fixed_agency$Pr[2], child_adu_test_fixed_pc$Pr[2],
           child_adu_test_fixed_wob$Pr[2], child_adu_test_fixed_gender$Pr[2], 
           child_adu_test_fixed_group$Pr[2], child_adu_test_EV_int$Pr[2], 
           child_adu_test_uncertainty_int$Pr[2], child_adu_test_agency_int$Pr[2])


no_compt_mods_fixed_effects_child_vs_adu.csv <- data.frame(variable, estimate, sem, chisq, df, p_value)

# Save fixed effects as a file
write.csv(no_compt_mods_fixed_effects_child_vs_adu.csv, sprintf("%s/no_compt_mods_fixed_effects_child_vs_adu.csv", wd), row.names = FALSE)

