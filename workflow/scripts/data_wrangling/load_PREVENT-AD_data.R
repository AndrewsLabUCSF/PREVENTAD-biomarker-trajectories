
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
diagnosis_file <- file.path(PREVENTAD_DATA_PATH,
                            "Clinical_diagnosis_Registered_PREVENTAD_2025-03-21.csv")
genetics_file <- file.path(PREVENTAD_DATA_PATH,
                           "Genetics_Registered_PREVENTAD_2025-03-21.csv")
lab_file <- file.path(PREVENTAD_DATA_PATH,
                      "Lab_Registered_PREVENTAD_2025-03-21.csv")
medical_history_file <- file.path(PREVENTAD_DATA_PATH,
                                  "EL_Medical_history_Registered_PREVENTAD_2025-03-21.csv")
plasma_4plex_file <- file.path(PREVENTAD_DATA_PATH,
                               "Plasma_4plex_IPMS_Registered_PREVENTAD_2025-03-21.csv")
plasma_ptau_file <- file.path(PREVENTAD_DATA_PATH,
                              "Plasma_p-tau217_Registered_PREVENTAD_2025-03-21.csv")
questionnaire_file <- file.path(PREVENTAD_DATA_PATH, 
                                "SelfReport_Behavioral_Questionnaires_Registered_PREVENTAD_2024-08-01.csv")

bp_pulse_weight <- read.csv(bp_pulse_weight_file)
demographics <- read.csv(demographics_file)
diagnosis <- read.csv(diagnosis_file) %>%
  select(CONP_ID, Clinical_diagnosis)
genetics <- read.csv(genetics_file)
lab <- read.csv(lab_file)
medical_history <- read.csv(medical_history_file)
plasma_4plex <- read.csv(plasma_4plex_file)
plasma_ptau <- read.csv(plasma_ptau_file)
questionnaire <- read.csv(questionnaire_file)


# COGDRISK --------------------------------------------------------------
## VARIABLE SELECTION ---------------------------------------------------

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
  left_join((lab %>%
               filter(!is.na(LDL_value)) %>%
               group_by(CONP_ID) %>%
               slice(1) %>%
               select(CONP_ID, LDL_value)),
            by="CONP_ID") %>%
  left_join((lab %>%
               select(CONP_ID, hba1c_value) %>%
               filter(!is.na(hba1c_value)) %>%
               group_by(CONP_ID) %>%
               slice(1)),
            by="CONP_ID") %>%
  left_join((medical_history %>%
               select(CONP_ID, past_depression, past_atrial_fibrillation, treatment_diabetes)),
            by="CONP_ID") %>%
  left_join((questionnaire %>%
               select(CONP_ID, head_injury_hospitalized, head_injury_severe) %>%
               filter_at(vars(head_injury_severe, head_injury_hospitalized), any_vars(!is.na(.)))),
            by="CONP_ID")


# Family history table
fhx_dataset <- demographics %>%
  select(CONP_ID, father_dx_ad_dementia:other_paternal_family_members_AD)


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


# APOE table
apoe_dataset <- genetics %>%
  select(CONP_ID, APOE)

# ptau table
ptau_dataset <- plasma_ptau %>%
  select(CONP_ID, Study_visit_label, Candidate_Age, ptau217_simoa_UGOT) %>%
  mutate(Candidate_Age = Candidate_Age/12)

# 4plex ipms table
ipms_dataset <- plasma_4plex %>%
  select(CONP_ID, Study_visit_label, Candidate_Age, AB_ratio_MS_UGOT, t_tau_simoa_UGOT,
         ptau181_simoa_UGOT, NFL_simoa_4plex, GFAP_simoa_4plex, AB_ratio_simoa_4plex) %>%
  mutate(Candidate_Age = Candidate_Age/12)


## RECODING ---------------------------------------------------------------

clinical_dataset_cogdrisk <- clinical_dataset %>%
  mutate(
    Sex = if_else(Sex == "Female", 0, 1) %>% factor(), # female == 0, male == 1
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
    High_cholesterol = if_else(total_cholesterol_value > 6.5, 1, 0) %>% factor(), # high == 1
    Depression = past_depression %>% factor(),
    Atrial_fibrillation = past_atrial_fibrillation %>% factor(),
    Diabetes = case_when(treatment_diabetes > 0 ~ 1,
                         Age < 60 & hba1c_value <= 0.061 ~ 0,
                         (Age >= 60 & Age < 70) & hba1c_value <= 0.075 ~ 0,
                         Age >= 70 & hba1c_value <= 0.07 ~ 0,
                         TRUE ~ 1) %>% factor(),
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
  relocate(Diabetes, .after=treatment_diabetes) %>%
  select(-c(Candidate_Age, Systolic_blood_pressure, Diastolic_blood_pressure,
            Height, Weight, total_cholesterol_value, past_depression,
            past_atrial_fibrillation, treatment_diabetes, hba1c_value,
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
         # social_engagement: 0 == not lonely, 1 == lonely
         Social_engagement = as.factor(if_else(Social_engagement_holder > mean(Social_engagement_holder, na.rm=TRUE), 
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
  # Physical inactivity: 0 = physically inactive, 1 = physically active
  mutate(Physical_inactivity = if_else((moderate_minutes_week + heavy_minutes_week) > 150, 0, 1)) %>%
  select(CONP_ID, Smoking, Cognitive_engagement, Social_engagement, Physical_inactivity)

# Merge
data_cogdrisk <- clinical_dataset_cogdrisk %>%
  left_join(lifestyle_dataset_cogdrisk, by="CONP_ID")



# LIBRA -------------------------------------------------------------------

libra_dataset <- clinical_dataset %>%
  select(CONP_ID, Systolic_blood_pressure, Diastolic_blood_pressure, total_cholesterol_value, 
         LDL_value, Height, Weight, hba1c_value, treatment_diabetes, past_depression) %>%
  left_join((lifestyle_dataset %>%
               select(CONP_ID, smoking_present)),
            by="CONP_ID") %>%
  
# Recoding
  mutate(
    Hypertension = if_else(
      Systolic_blood_pressure >= 140 | Diastolic_blood_pressure >= 90, 1, 0),
    Hypercholesterolemia = if_else(
      total_cholesterol_value >= 5.2 & LDL_value >= 4.2, 1, 0),
    BMI = Weight/(Height/100)^2,
    BMI_category = case_when(BMI < 18.5 ~ "Underweight",
                             BMI > 18.5 & BMI < 25 ~ "Normal",
                             BMI >= 25 & BMI < 30 ~ "Overweight",
                             BMI >= 30 ~ "Obese",
                             TRUE ~ NA) %>% factor(),
    Diabetes = case_when(treatment_diabetes > 0 ~ 1,
                         hba1c_value >= 0.065 ~ 1,
                         TRUE ~ 0),
    Depression = past_depression %>% factor(),
    Smoking = if_else(smoking_present >= 3, 1, 0)
    ) %>%
  relocate(Hypertension, .after=Diastolic_blood_pressure) %>%
  relocate(Hypercholesterolemia, .after=LDL_value) %>%
  relocate(BMI, .after=Weight) %>%
  relocate(BMI_category, .after=BMI) %>%
  relocate(Diabetes, .after=BMI_category) %>%
  relocate(Depression, .after=past_depression) %>%
  relocate(Smoking, .after=smoking_present) %>%
  select(-c(Systolic_blood_pressure, Diastolic_blood_pressure, total_cholesterol_value, 
            LDL_value, Height, Weight, hba1c_value, treatment_diabetes,
            past_depression, smoking_present))

