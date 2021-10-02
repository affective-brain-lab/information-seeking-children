# AUTHOR: GAIA MOLINARO
# Data analysis for supplementary experiment 3
# in Molinaro, Cogliati Dezza, & Sharot (in prep.)

# SET WORKING DIRECTORY
wd = "insert_your_wd_here"

#### SUPPLEMENTARY EXPERIMENT 3 (INSTRUMENTAL UTILITY CONTROL) ####
#### SE3 Load and prep datasets #### 
setwd(wd)
dat_ins_ctrl <- read.csv("lookit_ins_control_updated_trials.csv")
dat_ins_ctrl <- subset(dat_ins_ctrl, !(condition %in% c("catch_1", "catch_2")))
dat_ins_ctrl <- subset(dat_ins_ctrl, !(age_in_years == 13))
dat_ins_ctrl <- subset(dat_ins_ctrl, catch_trials_score == 100)

# Create age groups
dat_ins_ctrl$age_group <- "None"
dat_ins_ctrl$age_group[dat_ins_ctrl$age_in_years  %in% c(4, 5)] <- "4-5"
dat_ins_ctrl$age_group[dat_ins_ctrl$age_in_years  %in% c(6, 7)] <- "6-7"
dat_ins_ctrl$age_group[dat_ins_ctrl$age_in_years  %in% c(8, 9)] <- "8-9"
dat_ins_ctrl$age_group[dat_ins_ctrl$age_in_years  %in% c(10, 11, 12)] <- "10-12"
dat_ins_ctrl$age_group <- factor(dat_ins_ctrl$age_group, levels=c("4-5", "6-7", "8-9","10-12"))
dat_ins_ctrl$age_group_coded <- "None"
dat_ins_ctrl$age_group_coded[dat_ins_ctrl$age_in_years  %in% c(4, 5)] <- 0
dat_ins_ctrl$age_group_coded[dat_ins_ctrl$age_in_years  %in% c(6, 7)] <- 1
dat_ins_ctrl$age_group_coded[dat_ins_ctrl$age_in_years  %in% c(8, 9)] <- 2
dat_ins_ctrl$age_group_coded[dat_ins_ctrl$age_in_years  %in% c(10, 11, 12)] <- 3
dat_ins_ctrl$age_group_coded <- factor(dat_ins_ctrl$age_group_coded, levels=c("4-5", "6-7", "8-9","10-12"))


# Turn variables into appropriate data types
# Use scale() for standardized values
dat_ins_ctrl$non_z_delta_agency <- (as.numeric(dat_ins_ctrl$delta_agency))
dat_ins_ctrl$delta_agency <- scale(as.numeric(dat_ins_ctrl$delta_agency))
dat_ins_ctrl$age_in_years <- scale(as.numeric(dat_ins_ctrl$age_in_years))
dat_ins_ctrl$age_group_coded <- as.factor(dat_ins_ctrl$age_group_coded)
dat_ins_ctrl$RT_info_choice <- scale(as.numeric(dat_ins_ctrl$RT_info_choice))
dat_ins_ctrl$gorilla_ID <- as.factor(dat_ins_ctrl$gorilla_ID)
dat_ins_ctrl$gender_coded <- as.factor(dat_ins_ctrl$gender_coded)
contrasts(dat_ins_ctrl$gender_coded) <- contr.helmert(length(table(dat_ins_ctrl$gender_coded)))
dat_ins_ctrl$info_choice <- as.factor(dat_ins_ctrl$info_choice) # 0 = left, 1 = right
dat_ins_ctrl$percent_comprehension_non_Z <- as.numeric(dat_ins_ctrl$percent_comprehension)
dat_ins_ctrl$percent_comprehension <- scale(as.numeric(dat_ins_ctrl$percent_comprehension))
dat_ins_ctrl$chance <- 0.5

#### SE2 Comprehension scores ####
mean(subset(dat_ins_ctrl, condition == "instrumental_1")$percent_comprehension_non_Z)


#### SE3 Models by age group ####
# Children
child_ins_ctrl_mod_full <- glmer(info_choice ~ delta_agency + 
                                   percent_comprehension + gender_coded + 
                                   (delta_agency
                                    | gorilla_ID), 
                                 data = dat_ins_ctrl, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

sjPlot::tab_model(child_ins_ctrl_mod_full, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)


#### SE3 Betas plot #### 
# Create betas plot
ins_ctrl_instrumental_beta <- c(summary(child_ins_ctrl_mod_full)$coefficients[2])
ins_ctrl_instrumental_sem_lower <- c(summary(child_ins_ctrl_mod_full)$coefficients[2] - summary(child_ins_ctrl_mod_full)$coefficients[2,2])
ins_ctrl_instrumental_sem_upper <- c(summary(child_ins_ctrl_mod_full)$coefficients[2] + summary(child_ins_ctrl_mod_full)$coefficients[2,2])
instrumental_color <- "#4DAF4A"

ins_ctrl_motive_betas <- ggplot()+
  scale_y_continuous(breaks=seq(0,2,0.1), limits=c(0, 2), expand=c(0, 0)) + 
  # Children
  geom_errorbar(aes(0.3, ymin = ins_ctrl_instrumental_sem_lower[1], ymax = ins_ctrl_instrumental_sem_upper[1]), size = 0.8, width = .05, colour = instrumental_color) +
  geom_point(aes(0.3, ins_ctrl_instrumental_beta[1]), size = 4, color = instrumental_color)+
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size=20)) +
  labs(x="Utilities", y="Standardized beta predicting information-seeking")
ins_ctrl_motive_betas
