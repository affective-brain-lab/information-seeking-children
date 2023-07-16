# AUTHOR: GAIA MOLINARO
# Data analysis for experiment 5
# in Molinaro, Cogliati Dezza, Buehler, Moutsiana & Sharot (in prep.)

#### FUNCTIONS ####
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


# SET WORKING DIRECTORY
wd = "insert_your_wd_here"

#### LOAD REQUIRED PACKAGES #####
packages <- c("psych", "optimx", "ggplot2", "lme4", "dplyr")
pacman::p_load(packages, character.only = TRUE)


#### EXPERIMENT 5 (UNCERTAINTY TRACKING) ####
#### Load and prep datasets #### 
setwd(wd)
dat_unc_tr_ctrl <- read.csv("data_experiment_5.csv")

# Turn variables into appropriate data types
# Use scale() for standardized values
dat_unc_tr_ctrl$answerL <- (as.numeric(dat_unc_tr_ctrl$answerL))
dat_unc_tr_ctrl$answerR <- (as.numeric(dat_unc_tr_ctrl$answerR))
dat_unc_tr_ctrl$confidenceL <- (as.numeric(dat_unc_tr_ctrl$confidenceL))
dat_unc_tr_ctrl$confidenceR <- (as.numeric(dat_unc_tr_ctrl$confidenceR))
dat_unc_tr_ctrl$uncertainty_level_L <- (as.numeric(dat_unc_tr_ctrl$uncertainty_level_L))
dat_unc_tr_ctrl$uncertainty_level_R <- (as.numeric(dat_unc_tr_ctrl$uncertainty_level_R))
dat_unc_tr_ctrl$age_in_years_non_Z <- as.numeric(dat_unc_tr_ctrl$age_in_years)
dat_unc_tr_ctrl$age_in_years <- scale(as.numeric(dat_unc_tr_ctrl$age_in_years))
dat_unc_tr_ctrl$age_group <- "None"
dat_unc_tr_ctrl$age_group[dat_unc_tr_ctrl$age_in_years_non_Z  %in% c(4, 5)] <- "4-5"
dat_unc_tr_ctrl$age_group[dat_unc_tr_ctrl$age_in_years_non_Z  %in% c(6, 7)] <- "6-7"
dat_unc_tr_ctrl$age_group[dat_unc_tr_ctrl$age_in_years_non_Z  %in% c(8, 9)] <- "8-9"
dat_unc_tr_ctrl$age_group[dat_unc_tr_ctrl$age_in_years_non_Z  %in% c(10, 11, 12)] <- "10-12"
dat_unc_tr_ctrl$age_group <- factor(dat_unc_tr_ctrl$age_group, levels=c("4-5", "6-7", "8-9","10-12"))
dat_unc_tr_ctrl$age_group_coded <- as.factor(dat_unc_tr_ctrl$age_group_coded)
dat_unc_tr_ctrl$subject_ID <- as.factor(dat_unc_tr_ctrl$subject_ID)
dat_unc_tr_ctrl$gender_coded <- as.factor(dat_unc_tr_ctrl$gender_coded)
contrasts(dat_unc_tr_ctrl$gender_coded) <- contr.helmert(2)
dat_unc_tr_ctrl$percent_comprehension_non_Z <- as.numeric(dat_unc_tr_ctrl$percent_comprehension)
dat_unc_tr_ctrl$percent_comprehension <- scale(as.numeric(dat_unc_tr_ctrl$percent_comprehension))

dat_unc_tr_ctrl$delta_confidence_non_Z <- as.numeric(dat_unc_tr_ctrl$confidenceR - dat_unc_tr_ctrl$confidenceL)
dat_unc_tr_ctrl$delta_unc_non_Z <- as.numeric(dat_unc_tr_ctrl$uncertainty_level_R - dat_unc_tr_ctrl$uncertainty_level_L)

dat_unc_tr_ctrl$delta_confidence <- scale(as.numeric(dat_unc_tr_ctrl$confidenceR - dat_unc_tr_ctrl$confidenceL), center=FALSE)
dat_unc_tr_ctrl$delta_unc <- scale(as.numeric(dat_unc_tr_ctrl$uncertainty_level_R - dat_unc_tr_ctrl$uncertainty_level_L), center=FALSE)

#### Demographics ####
summary(subset(dat_unc_tr_ctrl, trial_number==1)$age_in_years_non_Z)
sdamr::sample_sd(subset(dat_unc_tr_ctrl, trial_number==1)$age_in_years_non_Z)
table(subset(dat_unc_tr_ctrl, trial_number==1)$gender)
table(subset(dat_unc_tr_ctrl, trial_number==1)$age_in_years_non_Z)

#### Comprehension scores ####
mean(subset(dat_unc_tr_ctrl, trial_number==1)$percent_comprehension_non_Z)
sdamr::sample_sd(subset(dat_unc_tr_ctrl, trial_number==1)$percent_comprehension_non_Z)
psych::describeBy(subset(dat_unc_tr_ctrl, trial_number==1)$percent_comprehension_non_Z, 
                  subset(dat_unc_tr_ctrl, trial_number==1)$age_in_years_non_Z)


#### LMEM ####
unc_tr_mod_full <- lmer(delta_confidence ~ delta_unc*age_in_years + percent_comprehension + (delta_unc|subject_ID),
                        data = dat_unc_tr_ctrl)
unc_tr_mod_drop_delta_unc_fixed <- lmer(delta_confidence ~ delta_unc:age_in_years + age_in_years + percent_comprehension + (delta_unc|subject_ID),
                        data = dat_unc_tr_ctrl)
unc_tr_mod_drop_delta_unc_int <- lmer(delta_confidence ~ delta_unc + age_in_years + percent_comprehension + (delta_unc|subject_ID),
                                        data = dat_unc_tr_ctrl)
anova(unc_tr_mod_full, unc_tr_mod_drop_delta_unc_fixed)
anova(unc_tr_mod_full, unc_tr_mod_drop_delta_unc_int)
sjPlot::tab_model(unc_tr_mod_full, show.se=TRUE)

#### Correlations by subject ####
corr_by_subj <- do.call( rbind, lapply( split(dat_unc_tr_ctrl, dat_unc_tr_ctrl$subject_ID),
                                        function(x) data.frame(subject_ID=x$subject_ID[1], corr=cor(x$delta_confidence_non_Z, x$delta_unc_non_Z, method="pearson")) ) )

dat_corr_by_subj <- data.frame(subject_ID=corr_by_subj$subject_ID, corr=corr_by_subj$corr)
dat_corr_by_subj <- merge(dat_corr_by_subj, 
                          select(subset(dat_unc_tr_ctrl, dat_unc_tr_ctrl$trial_number==1), 
                                 c("subject_ID", "age_in_years", "age_in_years_non_Z", "age_group_coded", "gender", "gender_coded", "percent_comprehension_non_Z")))

dat_corr_by_subj$min_corr = -dat_corr_by_subj$corr

psych::describeBy(dat_corr_by_subj$min_corr,  dat_corr_by_subj$age_in_years_non_Z)
aggregate(dat_corr_by_subj$min_corr, list(dat_corr_by_subj$age_in_years), FUN=mean, na.rm=TRUE)$x


#### Johnson-Neyman plot ####
n_tests <- (12*max(dat_unc_tr_ctrl$age_in_years_non_Z)) - (12*min(dat_unc_tr_ctrl$age_in_years_non_Z))
slopes_unc_tr <- modelbased::estimate_slopes(unc_tr_mod_full, trend = "delta_unc", at = "age_in_years", length = n_tests)

x_pos <- c(4, 5, 6, 7, 8, 9, 10, 11, 12)
x_pos_Z <- find_Z(x_pos, dat_unc_tr_ctrl$age_in_years_non_Z)
lims <- find_Z(c(4,12), dat_unc_tr_ctrl$age_in_years_non_Z)

# plot significance of the effect across time (Johnson-Neyman plot)
# modelbased package
jn_plot_unc_tr <- plot(slopes_unc_tr) + 
  #scale_fill_manual(values=c(nonsig_color, unc_color)) + 
  scale_fill_manual(values=c(unc_color)) + 
  theme_classic() +
  scale_y_continuous(limits = c(-1.5, 0.5)) +
  scale_x_continuous(breaks=x_pos_Z, labels=x_pos, limits=lims) +
  labs(x="Age (years)", y="Estimated β relating Δuncertainty\nto confidence", title="") + 
  theme(text = element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))
jn_plot_unc_tr
slopes_unc_tr

#### Load data from the main study #### 
#### Load and prep datasets #### 
rescale_variables <- TRUE  # whether to rescale the numeric variables
center_variables <- TRUE # whether to center numeric variables (except deltas)
rescale_age <- TRUE  # whether to rescale age
center_age <- TRUE # whether to center age

# Whether to use covariates in the model
use_covariates <- TRUE

###### Load  data ###### 
data_batch <- "" # either 1 or 2
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
dat$age_group_coded <- as.factor(dat$age_group_coded)
dat$subject_ID <- as.factor(dat$subject_ID)
dat$gender_coded[dat$gender_coded==3] <- 1 # 1 = male, 2 = female, 3 = other
dat$gender_coded <- as.factor(dat$gender_coded)
contrasts(dat$gender_coded) <- contr.helmert(2)
dat$info_choice <- as.factor(dat$info_choice) # 0 = left, 1 = right
dat$percent_comprehension_non_Z <- as.numeric(dat$percent_comprehension)
dat$percent_comprehension <- scale(as.numeric(dat$percent_comprehension), center=center_variables, scale=rescale_variables)
dat$wob_non_Z <- as.numeric(dat$wob)
dat$wob <- scale(as.numeric(dat$wob), center=center_variables, scale=rescale_variables)
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


#### Controlling for average uncertainty tracking in the main experiment ####

unc_tr <- c(aggregate(dat_corr_by_subj$min_corr, list(dat_corr_by_subj$age_in_years), FUN=mean, na.rm=TRUE)$x)
dat$unc_tr <- "None"
dat$unc_tr[dat$age_in_years_non_Z  == 4] <- unc_tr[1]
dat$unc_tr[dat$age_in_years_non_Z  == 5] <- unc_tr[2]
dat$unc_tr[dat$age_in_years_non_Z  == 6] <- unc_tr[3]
dat$unc_tr[dat$age_in_years_non_Z  == 7] <- unc_tr[4]
dat$unc_tr[dat$age_in_years_non_Z  == 8] <- unc_tr[5]
dat$unc_tr[dat$age_in_years_non_Z  == 9] <- unc_tr[6]
dat$unc_tr[dat$age_in_years_non_Z  == 10] <- unc_tr[7]
dat$unc_tr[dat$age_in_years_non_Z  == 11] <- unc_tr[8]
dat$unc_tr[dat$age_in_years_non_Z  == 12] <- unc_tr[9]
dat$unc_tr <- scale(as.numeric(dat$unc_tr))


{
  child_mod_with_unc_tr <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                           percent_comprehension + wob + fishing + unc_tr +
                           (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                         data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)

  child_mod_with_unc_tr_drop_Af_fixed <- glmer(info_choice ~ (delta_uncertainty_level + delta_agency)*age_in_years + delta_EV:age_in_years +
                                        percent_comprehension + wob + fishing + unc_tr +
                                        (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                      data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_with_unc_tr_drop_C_fixed <- glmer(info_choice ~ (delta_EV + delta_agency)*age_in_years + delta_uncertainty_level:age_in_years +
                                                       percent_comprehension + wob + fishing + unc_tr +
                                                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_with_unc_tr_drop_Ac_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level)*age_in_years + delta_agency:age_in_years +
                                                       percent_comprehension + wob + fishing + unc_tr +
                                                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_with_unc_tr_drop_Af_int <- glmer(info_choice ~ (delta_uncertainty_level + delta_agency)*age_in_years + delta_EV +
                                                       percent_comprehension + wob + fishing + unc_tr +
                                                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_with_unc_tr_drop_C_int <- glmer(info_choice ~ (delta_EV + delta_agency)*age_in_years + delta_uncertainty_level + 
                                                       percent_comprehension + wob + fishing + unc_tr +
                                                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_with_unc_tr_drop_Ac_int <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level)*age_in_years + delta_agency + 
                                                       percent_comprehension + wob + fishing + unc_tr +
                                                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_with_unc_tr_drop_unc_tr_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                                       percent_comprehension + wob + fishing +
                                                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_with_unc_tr_drop_age_in_years_fixed <- glmer(info_choice ~ delta_EV + delta_uncertainty_level + delta_agency + (delta_EV + delta_uncertainty_level + delta_agency):age_in_years +
                                                       percent_comprehension + wob + fishing + unc_tr +
                                                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_with_unc_tr_drop_pc_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                                        wob + fishing + unc_tr +
                                                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_with_unc_tr_drop_wob_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                                       percent_comprehension + fishing + unc_tr +
                                                       (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                                     data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
  child_mod_with_unc_tr_drop_fishing_fixed <- glmer(info_choice ~ (delta_EV + delta_uncertainty_level + delta_agency)*age_in_years +
                                        percent_comprehension + wob + unc_tr +
                                        (delta_EV + delta_uncertainty_level + delta_agency | subject_ID),
                                      data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=2e5)), nAGQ = 0)
  
}


test_child_with_unc_tr_EV_int <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_Af_int, test="Chisq")
test_child_with_unc_tr_uncertainty_int <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_C_int, test="Chisq")
test_child_with_unc_tr_agency_int <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_Ac_int, test="Chisq")
test_child_with_unc_tr_EV_fixed <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_Af_fixed, test="Chisq")
test_child_with_unc_tr_uncertainty_fixed <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_C_fixed, test="Chisq")
test_child_with_unc_tr_agency_fixed <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_Ac_fixed, test="Chisq")
test_child_with_unc_tr_age_in_years_fixed <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_age_in_years_fixed, test="Chisq")
test_child_with_unc_tr_pc_fixed <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_pc_fixed, test="Chisq")
test_child_with_unc_tr_wob_fixed <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_wob_fixed, test="Chisq")
test_child_with_unc_tr_fishing_fixed <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_fishing_fixed, test="Chisq")
test_child_with_unc_tr_unc_tr_fixed <- anova(child_mod_with_unc_tr, child_mod_with_unc_tr_drop_unc_tr_fixed, test="Chisq")


sjPlot::tab_model(child_mod_with_unc_tr, transform = NULL, auto.label = FALSE, show.stat = TRUE, show.ci=FALSE, show.se=TRUE)

