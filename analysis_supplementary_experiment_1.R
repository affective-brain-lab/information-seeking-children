# AUTHOR: GAIA MOLINARO
# Data analysis for the supplementary experiment 1
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

#### SUPPLEMENTARY EXPERIMENT 1 (NO COMPETITION) ####
#### SE1 Load and prep datasets #### 
setwd("C:/Gaia/ABL")
dat_no_compt <- read.csv("lookit_3_updated_trials.csv")
dat_no_compt <- subset(dat_no_compt, !(condition %in% c("catch_1", "catch_2")))
dat_no_compt <- subset(dat_no_compt, !(age_in_years == 13))
dat_no_compt <- subset(dat_no_compt, catch_trials_score == 100)
dat_no_compt <- subset(dat_no_compt, wob > 0.5)
# Rename variables with mislabelled heading
dat_no_compt$delta_uncertainty_level <- dat_no_compt$delta_uncertainty
dat_no_compt$gorilla_ID <- dat_no_compt$ï..gorilla_ID

# Create age groups
dat_no_compt$age_group <- "None"
dat_no_compt$age_group[dat_no_compt$age_in_years  %in% c(4, 5)] <- "4-5"
dat_no_compt$age_group[dat_no_compt$age_in_years  %in% c(6, 7)] <- "6-7"
dat_no_compt$age_group[dat_no_compt$age_in_years  %in% c(8, 9)] <- "8-9"
dat_no_compt$age_group[dat_no_compt$age_in_years  %in% c(10, 11, 12)] <- "10-12"
dat_no_compt$age_group <- factor(dat_no_compt$age_group, levels=c("4-5", "6-7", "8-9","10-12"))

# Turn variables into appropriate data types
# Use scale() for standardized values
dat_no_compt$non_z_delta_EV <- (as.numeric(dat_no_compt$delta_EV))
dat_no_compt$non_z_delta_uncertainty_level <- (as.numeric(dat_no_compt$delta_uncertainty_level))
dat_no_compt$non_z_delta_agency<- (as.numeric(dat_no_compt$delta_agency))
dat_no_compt$delta_EV <- scale(as.numeric(dat_no_compt$delta_EV))
dat_no_compt$delta_uncertainty_level <- scale(as.numeric(dat_no_compt$delta_uncertainty_level))
dat_no_compt$delta_agency <- scale(as.numeric(dat_no_compt$delta_agency))
dat_no_compt$age_in_years <- scale(as.numeric(dat_no_compt$age_in_years))
dat_no_compt$age_group_coded <- as.factor(dat_no_compt$age_group_coded)
dat_no_compt$RT_info_choice <- scale(as.numeric(dat_no_compt$RT_info_choice))
dat_no_compt$gorilla_ID <- as.factor(dat_no_compt$gorilla_ID)
dat_no_compt$gender_coded <- as.factor(dat_no_compt$gender_coded)
contrasts(dat_no_compt$gender_coded) <- contr.helmert(length(table(dat_no_compt$gender_coded)))
dat_no_compt$info_choice <- as.factor(dat_no_compt$info_choice) # 0 = left, 1 = right
dat_no_compt$percent_comprehension_non_Z <- as.numeric(dat_no_compt$percent_comprehension)
dat_no_compt$percent_comprehension <- scale(as.numeric(dat_no_compt$percent_comprehension))
dat_no_compt$wob_non_Z <- as.numeric(dat_no_compt$wob)
dat_no_compt$wob <- scale(as.numeric(dat_no_compt$wob))
dat_no_compt$chance <- 0.5

## SE1 - adults
setwd(wd)
adu_dat_no_compt <- read.csv("adults_4_updated_trials.csv")
adu_dat_no_compt <- subset(adu_dat_no_compt, !(condition %in% c("catch_1", "catch_2")))
adu_dat_no_compt <- subset(adu_dat_no_compt, (catch_trials_score == 100) & (wob > 0.5))

# Turn variables into appropriate data types
adu_dat_no_compt$delta_EV <- scale(as.numeric(adu_dat_no_compt$delta_EV))
adu_dat_no_compt$delta_uncertainty_level <- scale(as.numeric(adu_dat_no_compt$delta_uncertainty_level))
adu_dat_no_compt$delta_agency <- scale(as.numeric(adu_dat_no_compt$delta_agency))
adu_dat_no_compt$age_in_years_non_z <- as.numeric(adu_dat_no_compt$age_in_years)
adu_dat_no_compt$age_in_years <- scale(as.numeric(adu_dat_no_compt$age_in_years))
adu_dat_no_compt$age_group <- "Adult"
adu_dat_no_compt$age_group_coded <- as.factor(4)
adu_dat_no_compt$RT_info_choice <- scale(as.numeric(adu_dat_no_compt$RT_info_choice))
adu_dat_no_compt$gorilla_ID <- as.factor(adu_dat_no_compt$gorilla_ID)
adu_dat_no_compt$gender_coded <- as.factor(adu_dat_no_compt$gender_coded)
contrasts(adu_dat_no_compt$gender_coded) <- contr.helmert(2)
adu_dat_no_compt$info_choice <- as.factor(adu_dat_no_compt$info_choice) # 0 = left, 1 = right
adu_dat_no_compt$percent_comprehension_non_Z <- as.numeric(adu_dat_no_compt$percent_comprehension)
adu_dat_no_compt$percent_comprehension <- scale(as.numeric(adu_dat_no_compt$percent_comprehension))
adu_dat_no_compt$wob_non_Z <- as.numeric(adu_dat_no_compt$wob)
adu_dat_no_compt$wob <- scale(as.numeric(adu_dat_no_compt$wob))
adu_dat_no_compt$chance <- 0.5

## SE1 Merged data sets (children + adults) 
cols <- c("gorilla_ID", "info_choice", 
          "delta_EV", "delta_uncertainty_level", "delta_agency", 
          "gender_coded", "wob", "percent_comprehension", 
          "age_group", "age_group_coded")
dat_no_compt_child_adu <- rbind(dplyr::select(dat_no_compt, all_of(cols)), 
                                dplyr::select(adu_dat_no_compt, all_of(cols)))

# Create new variable and contrasts for children vs adults factor
dat_no_compt_child_adu$group <- as.factor(ifelse(dat_no_compt_child_adu$age_group_coded == 4, "Adults", "Children"))
contrasts(dat_no_compt_child_adu$group) <- contr.helmert(2)

#### SE1 Comprehension and wob scores ####
## t-tests

# Within group
t.test(subset(dat_no_compt, condition == "affective_1")$wob_non_Z, mu=0.5)
t.test(subset(adu_dat_no_compt, condition == "affective_1")$wob_non_Z, mu=0.5)

# Between groups
t.test(subset(dat_no_compt, condition == "affective_1")$percent_comprehension_non_Z, 
       subset(adu_dat_no_compt, condition == "affective_1")$percent_comprehension_non_Z)
t.test(subset(dat_no_compt, condition == "affective_1")$wob_non_Z, 
       subset(adu_dat_no_compt, condition == "affective_1")$wob_non_Z)

## Create long form data sets

# Create separate data sets with variables for children and adults' pc and wob
dat_no_compt_pc_wob <- dplyr::select(subset(dat_no_compt, condition == "affective_1"), 
                                     all_of(c("gorilla_ID", "age_group_coded", 
                                              "percent_comprehension_non_Z", 
                                              "wob_non_Z")))
dat_no_compt_pc_wob$percent_comprehension_child <- dat_no_compt_pc_wob$"percent_comprehension_non_Z"
dat_no_compt_pc_wob$wob_child <- dat_no_compt_pc_wob$"wob_non_Z"
adu_dat_no_compt_pc_wob <- dplyr::select(subset(adu_dat_no_compt,condition == "affective_1"), 
                                         all_of(c("gorilla_ID", "age_group_coded", 
                                                  "percent_comprehension_non_Z", 
                                                  "wob_non_Z")))

# Merge datasets
dat_no_compt_child_adu_pc_wob <- rbind(dplyr::select(dat_no_compt_pc_wob, all_of(
  c("gorilla_ID", "age_group_coded", "percent_comprehension_non_Z", "wob_non_Z"))), 
  dplyr::select(adu_dat_no_compt_pc_wob, all_of(
    c("gorilla_ID", "age_group_coded", "percent_comprehension_non_Z", "wob_non_Z"))))

# Create group factor 
dat_no_compt_child_adu_pc_wob$group <- as.factor(ifelse(dat_no_compt_child_adu_pc_wob$age_group_coded == 4, "Adults", "Children"))

# Turn into long form 
dat_no_compt_child_adu_pc_wob <- dat_no_compt_child_adu_pc_wob %>%
  tidyr::gather(key = "variable", value = "score", percent_comprehension_non_Z:wob_non_Z, factor_key=TRUE)

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
  theme_classic() + labs(x="comprehension", y="Score") + theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm"))
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
  theme_classic() + labs(x="comprehension", y="Score") + theme(legend.position="none", plot.margin=margin(0, 0, 0, 0, "cm"))
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
# Children
child_no_compt_mod_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                   percent_comprehension + wob + gender_coded + 
                                   (delta_EV + delta_uncertainty_level + delta_agency 
                                    | gorilla_ID), 
                                 data = dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)


# Adults
adu_no_compt_mod_full <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                 percent_comprehension + wob + gender_coded + 
                                 (delta_EV + delta_uncertainty_level + delta_agency 
                                  | gorilla_ID), 
                               data = adu_dat_no_compt, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)


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
  scale_y_continuous(breaks=seq(0,2.5,0.1), limits=c(0, 2.5), expand=c(0, 0)) + 
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
# Model
no_compt_mod_child_adu <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency +
                                  percent_comprehension + wob + gender_coded + group + 
                                  group:delta_EV + group:delta_uncertainty_level + group:delta_agency +
                                  (delta_EV + delta_uncertainty_level + delta_agency | gorilla_ID), 
                                data = dat_no_compt_child_adu, family = binomial, 
                                control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

# Model results
sjPlot::tab_model(no_compt_mod_child_adu, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)

