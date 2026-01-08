
# SETUP -------------------------------------------------------------------

# Load files
source('workflow/scripts/config.R')
source(CRS_FN)

clinical_raw <- readRDS(file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_clinical_imp_raw.rds"))
lifestyle_raw <- readRDS(file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_lifestyle_imp_raw.rds"))
meduse <- readRDS(file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_meduse.rds"))

## RECODING ---------------------------------------------------------------
clinical_cogd <- clinical_raw %>%
  left_join(meduse, by="CONP_ID") %>%
  left_join((lifestyle_raw %>% select(CONP_ID, gds_score, epworth_score, pittsburgh_total_score)),
            by="CONP_ID") %>%
  mutate(
    Sex = if_else(Sex == "Female", 0, 1) %>% factor(), # female == 0, male == 1
    Age = Candidate_Age/12,
    Education_level = case_when(
      Education_years > 11 ~ "High",
      Education_years > 8 & Education_years <= 11 ~ "Middle",
      Education_years < 8 ~ "Low",
      TRUE ~ NA) %>% factor(),
    
    # Hypertension risk factor present == 1
    Hypertension = case_when(
      medusage_antihypertensive == 1 ~ 1,
      Systolic_blood_pressure > 129 ~ 1,
      Diastolic_blood_pressure >= 80 ~ 1,
      medusage_antihypertensive == 0 & !is.na(Systolic_blood_pressure) & !is.na(Diastolic_blood_pressure) ~ 0,
      TRUE ~ NA) %>% factor(),
    
    BMI = Weight/(Height/100)^2,
    BMI_category = case_when(
      BMI < 18.5 ~ "Underweight",
      BMI > 18.5 & BMI < 25 ~ "Normal",
      BMI >= 25 & BMI < 30 ~ "Overweight",
      BMI >= 30 ~ "Obese",
      TRUE ~ NA) %>% factor(),
    
    # High cholesterol risk factor present == 1
    High_cholesterol = if_else(total_cholesterol_value > 6.5, 1, 0) %>% factor(), 
    
    # Insomnia risk factor present == 1
    Insomnia = if_else(
      epworth_score >= 10 & pittsburgh_total_score > 5, 1, 0
      ) %>% factor(),
    
    # Depression risk factor present == 1
    Depression = case_when(
      past_depression == 1 ~ 1,
      medusage_antidepressants == 1 ~ 1,
      gds_score >= 5 ~ 1,
      past_depression == 0 & medusage_antidepressants == 0 ~ 0,
      TRUE ~ NA) %>% factor(),
    Atrial_fibrillation = past_atrial_fibrillation %>% factor(),
    
    # Diabetes risk factor present == 1
    Diabetes = case_when(
      treatment_diabetes > 0 ~ 1,
      medusage_antidiabetic == 1 ~ 1,
      Age < 60 & treatment_diabetes == 0 & hba1c_value <= 0.061 ~ 0,
      (Age >= 60 & Age < 70) & hba1c_value <= 0.075 & treatment_diabetes == 0 ~ 0,
      Age >= 70 & treatment_diabetes == 0 & hba1c_value <= 0.07 ~ 0,
      is.na(treatment_diabetes) & is.na(medusage_antidiabetic) & is.na(hba1c_value) ~ NA,
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
            past_atrial_fibrillation, hba1c_value, head_injury_hospitalized, 
            head_injury_severe, starts_with("medusage_"), starts_with("treatment_"),
            ends_with("medication"), diagnosed_impairment, epworth_score,
            gds_score, pittsburgh_total_score))


lifestyle_cogd <- lifestyle_raw %>%
  # smoking: 2 == current, 1 == former, 0 == non smoker
  mutate(Smoking = case_when(smoking_present == 4 ~ 2,
                             smoking_present == 3 ~ 2,
                             smoking_present == 2 ~ 1,
                             smoking_present == 1 ~ 1,
                             smoking_present == 0 ~ 0,
                             TRUE ~ NA),
         Cognitive_engagement = ntile(epoch_score_currently, 3) %>% factor(),
         Social_engagement_holder = rowMeans(select(., social_life_frequency_activities:social_life_frequency_phone_calls)),
         # social_engagement: 0 == not lonely, 1 == lonely
         Social_engagement = as.factor(if_else(Social_engagement_holder > mean(Social_engagement_holder, na.rm=TRUE), 
                                               0, 1))) %>%
  select(-Social_engagement_holder)

# Recoding exercise columns
lifestyle_exercise_cogd <- lifestyle_cogd %>%
  mutate(
    light_minutes_week = 0,
    moderate_minutes_week = 0,
    heavy_minutes_week = 0) %>%
  mutate(across(contains("_intensity"), ~ as.numeric(as.character(.))),
         across(contains("_days"), ~ as.numeric(as.character(.))),
         across(contains("_hours"), ~ as.numeric(as.character(.)))) %>%
  rowwise() %>%
  mutate(
    has_any_exercise_data = any(!is.na(c_across(contains("exer_curr_act"))))
  ) %>%
  ungroup()

for(i in 1:5) {  
  intensity_col <- paste0("exer_curr_act", i, "_intensity")
  days_col <- paste0("exer_curr_act", i, "_days")
  hours_col <- paste0("exer_curr_act", i, "_hours")
  
  lifestyle_exercise_cogd <- lifestyle_exercise_cogd %>%
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

lifestyle_exercise_cogd <- lifestyle_exercise_cogd %>%
  mutate(
    light_minutes_week = if_else(has_any_exercise_data, light_minutes_week, NA_real_),
    moderate_minutes_week = if_else(has_any_exercise_data, moderate_minutes_week, NA_real_),
    heavy_minutes_week = if_else(has_any_exercise_data, heavy_minutes_week, NA_real_)
  ) %>%
  select(-has_any_exercise_data)

lifestyle_cogd <- lifestyle_exercise_cogd %>%
  relocate(Smoking, .after=smoking_present) %>%
  relocate(Cognitive_engagement, .after=epoch_score_currently) %>%
  relocate(light_minutes_week, .after=Cognitive_engagement) %>%
  relocate(moderate_minutes_week, .after=light_minutes_week) %>%
  relocate(heavy_minutes_week, .after=moderate_minutes_week) %>%
  # Physical inactivity: 0 = physically active, 1 = physically inactive
  mutate(Physical_inactivity = if_else((moderate_minutes_week + heavy_minutes_week) > 150, 0, 1)) %>%
  select(CONP_ID, Smoking, Cognitive_engagement, Social_engagement, Physical_inactivity)

# Merge
dat_cogd <- clinical_cogd %>%
  left_join(lifestyle_cogd, by="CONP_ID")

# Scoring
dat_cogd_scored <- dat_cogd %>%
  mutate(score_cogd = calculate_cogdrisk(.)) %>%
  relocate(score_cogd, .after=CONP_ID) 

# Save
saveRDS(dat_cogd, 
        file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_cogd_dat.rds"))
saveRDS(dat_cogd_scored, 
        file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_cogd_scored_dat.rds"))
