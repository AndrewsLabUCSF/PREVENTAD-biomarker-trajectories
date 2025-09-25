
# SETUP -------------------------------------------------------------------

# Load configuration file
source('workflow/scripts/config.R')

# Path to data
PREVENTAD_DATA_PATH <- 'data/PREVENT-AD'

# Load datasets
bp_pulse_weight_file <- file.path(PREVENTAD_DATA_PATH, 
                                  "BP_Pulse_Weight_Registered_PREVENTAD_2025-03-21.csv")
demographics_file <- file.path(PREVENTAD_DATA_PATH, 
                               "Demographics_Registered_PREVENTAD_2025-03-21.csv")
medical_history_file <- file.path(PREVENTAD_DATA_PATH,
                                  "EL_Medical_history_Registered_PREVENTAD_2025-03-21.csv")
lab_file <- file.path(PREVENTAD_DATA_PATH,
                      "Lab_Registered_PREVENTAD_2025-03-21.csv")
questionnaire_file <- file.path(PREVENTAD_DATA_PATH, 
                                "SelfReport_Behavioral_Questionnaires_Registered_PREVENTAD_2024-08-01.csv")

bp_pulse_weight <- read.csv(bp_pulse_weight_file)
demographics <- read.csv(demographics_file)
lab <- read.csv(lab_file)
medical_history <- read.csv(medical_history_file)
questionnaire <- read.csv(questionnaire_file)


# VARIABLE SELECTION ------------------------------------------------------

# Clinical table
clinical_dataset <- demographics %>%
  select(CONP_ID, Sex, Education_years, Height) %>%
  left_join((bp_pulse_weight %>% 
               filter(str_detect(Visit_label, "^BL")) %>%
               select(CONP_ID, Candidate_Age, Systolic_blood_pressure, Diastolic_blood_pressure)), 
            by="CONP_ID") %>%
  left_join((bp_pulse_weight %>%
               group_by(CONP_ID) %>%
               filter(!is.na(Weight)) %>%
               slice(1) %>%
               select(CONP_ID, Weight)),
            by="CONP_ID") %>%
  left_join((lab %>%
               filter(!is.na(total_cholesterol_value)) %>%
               group_by(CONP_ID) %>%
               slice(1) %>%
               select(CONP_ID, total_cholesterol_value)),
            by="CONP_ID") %>%
  left_join((medical_history %>%
               select(CONP_ID, past_depression, past_atrial_fibrillation, treatment_diabetes)),
            by="CONP_ID") %>%
  left_join((questionnaire %>%
               select(CONP_ID, head_injury_hospitalized, head_injury_severe) %>%
               filter_at(vars(head_injury_severe, head_injury_hospitalized), any_vars(!is.na(.)))),
            by="CONP_ID")


# Lifestyle table
lifestyle_dataset <- demographics %>%
  select(CONP_ID) %>%
  left_join((questionnaire %>%
               select(CONP_ID, smoking_present) %>%
               filter(!is.na(smoking_present)) %>%
               group_by(CONP_ID) %>%
               slice_tail()),
            by="CONP_ID") %>%
  left_join(questionnaire %>%
              select(CONP_ID, epoch_score_currently) %>%
              filter(!is.na(epoch_score_currently)) %>%
              group_by(CONP_ID) %>%
              slice(1),
            by="CONP_ID") %>%
  left_join((questionnaire %>%
               select(CONP_ID, starts_with("exer_curr_")) %>%
               select(!matches("category")) %>%
               filter(!is.na(exer_curr_act1_intensity)) %>%
               group_by(CONP_ID) %>%
               slice(1)),
            by="CONP_ID") %>%
  left_join((questionnaire %>%
               select(CONP_ID, starts_with("social_life_")) %>%
               group_by(CONP_ID) %>%
               filter_at(vars(starts_with("social_life_")), any_vars(!is.na(.))) %>%
               slice(1)),
            by="CONP_ID")


# RECODING ----------------------------------------------------------------

clinical_dataset_cogdrisk <- clinical_dataset %>%
  mutate(
    Age = Candidate_Age/12,
    Education_level = case_when(Education_years > 11 ~ "High",
                                Education_years > 8 & Education_years <= 11 ~ "Middle",
                                Education_years < 8 ~ "Low",
                                TRUE ~ NA) %>% factor(),
    Hypertension = if_else(Systolic_blood_pressure > 129 | Diastolic_blood_pressure >= 80,
                           1, 0) %>% factor(),
    BMI = Weight/(Height/100)^2,
    BMI_category = case_when(BMI < 18.5 ~ "Underweight",
                             BMI > 18.5 & BMI < 25 ~ "Normal",
                             BMI >= 25 & BMI < 30 ~ "Overweight",
                             BMI >= 30 ~ "Obese",
                             TRUE ~ NA) %>% factor(),
    High_cholesterol = if_else(total_cholesterol_value > 6.5, 1, 0) %>% factor(),
    Depression = past_depression %>% factor(),
    Atrial_fibrillation = past_atrial_fibrillation %>% factor(),
    Diabetes_treatment = case_when(treatment_diabetes == 0 ~ "Never treated",
                                   treatment_diabetes == 1 ~ "Treated in the past",
                                   treatment_diabetes == 2 ~ "Currently treated",
                                   TRUE ~ NA) %>% factor(),
    TBI = if_else(head_injury_hospitalized == 1 | head_injury_severe == 1,
                  1, 0) %>% factor()) %>%
  relocate(Age, .after=Sex) %>%
  relocate(Education_level, .after=Education_years) %>%
  relocate(BMI, .after=Weight) %>%
  relocate(BMI_category, .after=BMI) %>%
  relocate(Hypertension, .after=BMI_category) %>%
  relocate(High_cholesterol, .after=total_cholesterol_value) %>%
  relocate(Depression, .after=past_depression) %>%
  relocate(Atrial_fibrillation, .after=past_atrial_fibrillation) %>%
  relocate(Diabetes_treatment, .after=treatment_diabetes) %>%
  select(-c(Candidate_Age, Systolic_blood_pressure, Diastolic_blood_pressure,
            Height, Weight, total_cholesterol_value, past_depression,
            past_atrial_fibrillation, treatment_diabetes,
            head_injury_hospitalized, head_injury_severe))


lifestyle_dataset_cogdrisk <- lifestyle_dataset %>%
  mutate(Smoking = case_when(smoking_present == 4 ~ 2,
                             smoking_present == 3 ~ 2,
                             smoking_present == 2 ~ 1,
                             smoking_present == 1 ~ 1,
                             smoking_present == 0 ~ 0,
                             TRUE ~ NA),
         Cognitive_engagement = ntile(epoch_score_currently, 3),
         Social_engagement_holder = rowMeans(select(., social_life_frequency_activities:social_life_frequency_phone_calls)),
         Social_engagement = as.factor(if_else(
           Social_engagement_holder > mean(lifestyle_dataset_cogdrisk$Social_engagement_holder, na.rm=TRUE), 
           0, 1))) %>%
  select(-Social_engagement_holder)
  
# Recoding exercise columns
lifestyle_exercise_cogdrisk <- lifestyle_dataset_cogdrisk %>%
  mutate(
    light_minutes_week = 0,
    moderate_minutes_week = 0,
    heavy_minutes_week = 0) %>%
  mutate(across(contains("_intensity"), ~ as.numeric(as.character(.))),
         across(contains("_days"), ~ as.numeric(as.character(.))),
         across(contains("_hours"), ~ as.numeric(as.character(.))))

for(i in 1:5) {  
  intensity_col <- paste0("exer_curr_act", i, "_intensity")
  days_col <- paste0("exer_curr_act", i, "_days")
  hours_col <- paste0("exer_curr_act", i, "_hours")
  
  lifestyle_exercise_cogdrisk <- lifestyle_exercise_cogdrisk %>%
    mutate(
      # Only calculate if ALL three values for this activity are present
      activity_minutes = case_when(
        !is.na(.data[[intensity_col]]) & 
          !is.na(.data[[days_col]]) & 
          !is.na(.data[[hours_col]]) ~ .data[[days_col]] * .data[[hours_col]] * 60,
        TRUE ~ 0),
      
      light_minutes_week = light_minutes_week + 
        ifelse(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 1, activity_minutes, 0),
      
      moderate_minutes_week = moderate_minutes_week + 
        ifelse(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 2, activity_minutes, 0),
      
      heavy_minutes_week = heavy_minutes_week + 
        ifelse(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 3, activity_minutes, 0)
    ) %>%
    select(-activity_minutes)
}

lifestyle_dataset_cogdrisk <- lifestyle_exercise_cogdrisk %>%
  relocate(Smoking, .after=smoking_present) %>%
  relocate(Cognitive_engagement, .after=epoch_score_currently) %>%
  relocate(light_minutes_week, .after=Cognitive_engagement) %>%
  relocate(moderate_minutes_week, .after=light_minutes_week) %>%
  relocate(heavy_minutes_week, .after=moderate_minutes_week) %>%
  select(-c(smoking_present)) %>%
  select(-matches("exer_curr_act"))

