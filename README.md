# information-seeking-children
Data and code from manuscript: Molinaro, Cogliati Dezza, Buehler, Moutsiana &amp; Sharot (in prep.)

## Data files
### data_experiment_1.csv (initial study) and data_experiment_2.csv (replication)
Contains the following variables:
* subject_ID: the participant's ID
* age_in_years: the participant's age in years
* age_in_months: the participant's age in months	
* gender: the participant's gender (Female, Male, Other)
* gender_coded: the participant's gender (Female=0, Male=1, Other=2)
* age_group: the participant's age group (4-5, 6-7, 8-9, 10-12)	
* age_group_coded: the participant's age group (4-5=0, 6-7=1, 8-9=2, 10-12=3)	
* info_choice: whether the participant wanted more information on the left or on the right side of the screen (left=0, right=1)	
* trial_number: the trial number
* condition: trial number (trial_n for normal trials, catch_n for catch trials)	
* catch_trials: choice on catch trials ("only cans"=0, "only fish"=1)	
* fishing_choiceL: whether the participant chose to fish on the left side, if given the option (yes=0, no=1)	
* fishing_choiceR: whether the participant chose to fish on the right side, if given the option (yes=0, no=1)	
* RT_catch: reaction time for catch_trials
* RT_info_choice: reaction time for info_choice	
* RT_fishingR: reaction time for fishing_choiceR
* RT_fishingL: reaction time for fishing_choiceR	
* items_countL: a dictionary containing each item's count for the left side	
* items_countR: a dictionary containing each item's count for the right side	 	
* agency_probL: the probability that the child would get to choose whether to fish on the left side (0, 50, 100)	
* agency_probR: the probability that the child would get to choose whether to fish on the right side (0, 50, 100)			
* EV_L: expected value on the left side	
* EV_R: expected value on the right side	
* SD_L: standard deviation on the left side	
* SD_R: standard deviation on the right side	
* uncertainty_level_L: number of seaweeds on the left side	
* uncertainty_level_R: number of seaweeds on the right side	
* delta_EV: difference in expected value between the two sides (right - left)	
* delta_SD: difference in standard deviation between the two sides (right - left)			
* delta_uncertainty_level: difference in number of seaweeds between the two sides	(right - left)		
* delta_agency: difference in agency probability between the two sides	(right - left)			
* items_count_flashL: a dictionary containing each item's count after using the flashlight on the left side	 	
* items_count_flashR: a dictionary containing each item's count after using the flashlight on the right side			
* wob: score on the EV comparison ("which one is better?") task	
* avg_wob_RT: average reaction time in the EV comparison ("which one is better?") task		
* percent_comprehension: score on the instructions comprehension task	
* q1-13: score on the question q_n in the instructions comprehension task (one column per question; 2=responded correctly immediately, 1=responded correctly after one try, 0=responded correctly after two tries)  	
* catch_trials_score: total score in the catch trials (2/2 correct=100, 1/2 correct=50, 0/2 corect=0)	
* trial_number: trial number	
* prop_correct_fishing: proportion of correct fishing choices, based on delta_EV
* fishing_rewardL: reward obtained after the fishing decision was made on the left
* fishing_rewardR: reward obtained after the fishing decision was made on the right

### data_experiment_3.csv
Contains the following variables:
* subject_ID: the participant's ID
* age_in_years: the participant's age in years
* age_in_months: the participant's age in months	
* gender: the participant's gender (Female, Male, Other)
* gender_coded: the participant's gender (Female=0, Male=1, Other=2)
* catch_trials_score: total score in the catch trials (2/2 correct=100, 1/2 correct=50, 0/2 corect=0)	
* percent_comprehension: score on the instructions comprehension task	
* info_choice: whether the participant wanted more information on the left or on the right side of the screen (left=0, right=1)	
* combination_order: whether side 1 corresponds to left and 2 to right (0) or vice-versa (1)	
* condition: condition (cognitive_n for normal trials, catch-n for catch trials)	
* catch_trials: choice on catch trials ("only cans"=0, "only fish"=1)	
* fishing_choiceL: whether the participant chose to fish on the left side, if given the option (yes=0, no=1)	
* fishing_choiceR: whether the participant chose to fish on the right side, if given the option (yes=0, no=1)	
* RT_catch: reaction time for catch_trials
* RT_info_choice: reaction time for info_choice	
* RT_fishingR: reaction time for fishing_choiceR
* RT_fishingL: reaction time for fishing_choiceR		
* agency_1: number of seaweeds on the first side (according to the predifined condition)	
* agency_2: number of seaweeds on the second side (according to the predifined condition)		
* agency_probL: number of seaweeds on the left side 	
* agency_probR: number of seaweeds on the right side 		
* delta_agency: difference in agency probability between the two sides (right - left)

The file is ordered by condition, although conditions appeared in random order in the experiment.

### data_experiment_4
Contains the following variables:
* subject_ID: the participant's ID
* age_in_years: the participant's age in years
* age_in_months: the participant's age in months	
* gender: the participant's gender (Female, Male, Other)
* gender_coded: the participant's gender (Female=0, Male=1, Other=2)
* catch_trials_score: total score in the catch trials (2/2 correct=100, 1/2 correct=50, 0/2 corect=0)	
* percent_comprehension: score on the instructions comprehension task	
* info_choice: whether the participant wanted more information on the left or on the right side of the screen (left=0, right=1)	
* combination_order: whether side 1 corresponds to left and 2 to right (0) or vice-versa (1)	
* condition: condition (cognitive_n for normal trials, catch-n for catch trials)	
* catch_trials: choice on catch trials ("only cans"=0, "only fish"=1)	
* fishing_choiceL: whether the participant chose to fish on the left side, if given the option (yes=0, no=1)	
* fishing_choiceR: whether the participant chose to fish on the right side, if given the option (yes=0, no=1)	
* RT_catch: reaction time for catch_trials
* RT_info_choice: reaction time for info_choice	
* RT_fishingR: reaction time for fishing_choiceR
* RT_fishingL: reaction time for fishing_choiceR		
* uncertainty_level_1: number of seaweeds on the first side (according to the predifined condition)	
* uncertainty_level_2: number of seaweeds on the second side (according to the predifined condition)		
* uncertainty_level_L: number of seaweeds on the left side 	
* uncertainty_level_R: number of seaweeds on the right side 		
* delta_uncertainty_level: difference in number of seaweeds between the two sides (right - left)

The file is ordered by condition, although conditions appeared in random order in the experiment.

### data_experiment_5
Contains the following variables:
* subject_ID: the participant's ID
* age_in_years: the participant's age in years
* age_in_months: the participant's age in months	
* gender: the participant's gender (Female, Male, Other)
* gender_coded: the participant's gender (Female=0, Male=1, Other=2)
* age_group: the participant's age group (4-5, 6-7, 8-9, 10-12)	
* age_group_coded: the participant's age group (4-5=0, 6-7=1, 8-9=2, 10-12=3)	
* trial_number: the trial number
* answerL: whether the participant thought the left-side fisherman would get a fish (0) or a can (1)	
* answerR: whether the participant thought the right-side fisherman would get a fish (0) or a can (1)	
* confidenceL: the participant's confidence level for the left-side decision
* confidenceR: the participant's confidence level for the left-side decision
* condition: trial number (trial_n for normal trials, catch_n for catch trials)	
* combination: whether the first fisherman's items appeared on the left (0) or the righ
* RT_answerL: reaction time for the left-side answer
* RT_answerR: reaction time for the right-side answer
* RT_confidenceL: reaction time for the left-side confidence
* RT_confidenceR: reaction time for the right-side confidence
* items_countL: a dictionary containing each item's count for the left side	
* items_countR: a dictionary containing each item's count for the right side	 		
* uncertainty_level_L: number of seaweeds on the left side	
* uncertainty_level_R: number of seaweeds on the right side	
* n_fish_L: number of fish on the left
* n_fish_R: number of fish on the right
* n_cans_L: number of cans on the left
* n_cans_R: number of cans on the right
* q1-5: score on the question q_n in the instructions comprehension task (one column per question; 2=responded correctly immediately, 1=responded correctly after one try, 0=responded correctly after two tries)  
* percent_comprehension: score on the instructions comprehension task	

### data_experiment_6.csv (adults initial study) and data_experiment_7.csv (adults replication)
Contains the same variables as data_main_experiment_1.csv, but without age_group and age_group_coded

## Data analysis files
### analysis_experiments_1_2.R
Contains the analyses for Experiment 1 and Experiment 2

### analysis_experiment_3.R
Contains the analyses for Experiment 3 

### analysis_experiment_4.R
Contains the analyses for Experiment 4

### analysis_experiment_5.R
Contains the analyses for Experiment 5

### analysis_experiments_6_7.R
Contains the analyses for Experiments 6 and 7
