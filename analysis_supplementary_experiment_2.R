# AUTHOR: GAIA MOLINARO
# Data analysis for the main and supplementary experiments
# in Molinaro, Cogliati Dezza, & Sharot (in prep.)

# SET WORKING DIRECTORY
# wd = "insert_your_wd_here"

#### SUPPLEMENTARY EXPERIMENT 2 (COGNITIVE UTILITY CONTROL) ####
#### SE2 Load and prep datasets #### 
setwd(wd)
dat_cog_ctrl <- read.csv("lookit_cog_control_updated_trials.csv")
dat_cog_ctrl <- subset(dat_cog_ctrl, !(condition %in% c("catch_1", "catch_2")))
dat_cog_ctrl <- subset(dat_cog_ctrl, !(age_in_years == 13))
dat_cog_ctrl <- subset(dat_cog_ctrl, catch_trials_score == 100)

# Rename variables with mislabelled heading
dat_cog_ctrl$delta_uncertainty_level <- dat_cog_ctrl$delta_uncertainty

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
dat_cog_ctrl$non_z_delta_uncertainty_level <- (as.numeric(dat_cog_ctrl$delta_uncertainty_level))
dat_cog_ctrl$delta_uncertainty_level <- scale(as.numeric(dat_cog_ctrl$delta_uncertainty_level))
dat_cog_ctrl$age_in_years <- scale(as.numeric(dat_cog_ctrl$age_in_years))
dat_cog_ctrl$age_group_coded <- as.factor(dat_cog_ctrl$age_group_coded)
dat_cog_ctrl$RT_info_choice <- scale(as.numeric(dat_cog_ctrl$RT_info_choice))
dat_cog_ctrl$gorilla_ID <- as.factor(dat_cog_ctrl$gorilla_ID)
dat_cog_ctrl$gender_coded <- as.factor(dat_cog_ctrl$gender_coded)
contrasts(dat_cog_ctrl$gender_coded) <- contr.helmert(length(table(dat_cog_ctrl$gender_coded)))
dat_cog_ctrl$info_choice <- as.factor(dat_cog_ctrl$info_choice) # 0 = left, 1 = right
dat_cog_ctrl$percent_comprehension_non_Z <- as.numeric(dat_cog_ctrl$percent_comprehension)
dat_cog_ctrl$percent_comprehension <- scale(as.numeric(dat_cog_ctrl$percent_comprehension))
dat_cog_ctrl$chance <- 0.5

#### SE2 Comprehension scores ####
mean(subset(dat_cog_ctrl, condition == "cognitive_1")$percent_comprehension_non_Z)

#### SE2 Models by age group ####
# Children
child_cog_ctrl_mod_full <- glmer(info_choice ~ delta_uncertainty_level + 
                                   percent_comprehension + gender_coded + 
                                   (delta_uncertainty_level
                                    | gorilla_ID), 
                                 data = dat_cog_ctrl, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

sjPlot::tab_model(child_cog_ctrl_mod_full, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)


#### SE2 Beta plot #### 
# Create beta plot
cog_ctrl_cognitive_beta <- c(summary(child_cog_ctrl_mod_full)$coefficients[2])
cog_ctrl_cognitive_sem_lower <- c(summary(child_cog_ctrl_mod_full)$coefficients[2] - summary(child_cog_ctrl_mod_full)$coefficients[2,2])
cog_ctrl_cognitive_sem_upper <- c(summary(child_cog_ctrl_mod_full)$coefficients[2] + summary(child_cog_ctrl_mod_full)$coefficients[2,2])
cognitive_color <- "#377EB8"

cog_ctrl_motive_betas <- ggplot()+
  scale_y_continuous(breaks=seq(0,2,0.1), limits=c(0, 2), expand=c(0, 0)) + 
  # Children
  geom_errorbar(aes(0.3, ymin = cog_ctrl_cognitive_sem_lower[1], ymax = cog_ctrl_cognitive_sem_upper[1]), size = 0.8, width = .05, colour = cognitive_color) +
  geom_point(aes(0.3, cog_ctrl_cognitive_beta[1]), size = 4, color = cognitive_color)+
 theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size=20)) +
  labs(x="Utilities", y="Standardized beta predicting information-seeking")
cog_ctrl_motive_betas
