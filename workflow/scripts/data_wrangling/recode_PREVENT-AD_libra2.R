
# SETUP -------------------------------------------------------------------

# Load files
source('workflow/scripts/config.R')

clinical_raw <- readRDS(file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_clinical_imp_raw.rds"))
lifestyle_raw <- readRDS(file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_lifestyle_imp_raw.rds"))


# RECODING ----------------------------------------------------------------
dat_libra <- clinical_raw %>%
  select(CONP_ID, Systolic_blood_pressure, Diastolic_blood_pressure, total_cholesterol_value, 
         LDL_value, Height, Weight, hba1c_value, treatment_diabetes, treatment_hypertension,
         treatment_hyperlipidemia, past_depression, diagnosed_impairment) %>%
  left_join((lifestyle_raw %>%
               select(CONP_ID, smoking_present, epoch_score_currently, pittsburgh_total_score,
                      social_life_frequency_activities:social_life_frequency_phone_calls,
                      gds_score)),
            by="CONP_ID") %>%
  mutate(
    # Sleep disturbances:
    # 1 == highest tertile
    Sleep_disturbances_holder = ntile(pittsburgh_total_score, 3),
    Sleep_disturbances = if_else(Sleep_disturbances_holder == 3, 1, 0),
    
    Hypertension = if_else(
      Systolic_blood_pressure >= 140 | Diastolic_blood_pressure >= 90 | treatment_hypertension == 2, 1, 0),
    
    # Hypercholesterolemia
    # 1 == TC >= 6.2 or currently treating hyperlipidemia
    Hypercholesterolemia = if_else(
      total_cholesterol_value >= 6.2 | treatment_hyperlipidemia == 2, 1, 0),
    BMI = Weight/(Height/100)^2,
    Obesity = if_else(BMI >= 30, 1, 0),
    
    # Hearing impairment:
    # 1 == yes 
    Hearing_impairment = if_else(diagnosed_impairment != 0, 1, 0),
    Diabetes = case_when(treatment_diabetes > 0 ~ 1,
                         hba1c_value >= 0.065 ~ 1,
                         TRUE ~ 0),
    
    # Depression
    # 1 == yes
    Depression = if_else(past_depression == 1 | gds_score >= 5, 1, 0),
    Smoking = if_else(as.numeric(smoking_present) >= 3, 1, 0),
    Cognitive_activity_holder = ntile(epoch_score_currently, 3),
    
    # Cognitive_activity:
    # 1 == highest tertile
    Cognitive_activity = if_else(Cognitive_activity_holder == 3, 1, 0)
  ) %>%
  
  # Social participation
  # 1 == lowest tertile
  mutate(Social_activity_total = rowSums(select(., social_life_frequency_activities:social_life_frequency_phone_calls)),
         Social_activity_holder = ntile(Social_activity_total, 3),
         Social_activity = if_else(Social_activity_holder == 1, 1, 0)) %>%
  
  relocate(Sleep_disturbances, .after=CONP_ID) %>%
  relocate(Hypertension, .after=Diastolic_blood_pressure) %>%
  relocate(Hypercholesterolemia, .after=LDL_value) %>%
  relocate(BMI, .after=Hypercholesterolemia) %>%
  relocate(Diabetes, .after=BMI) %>%
  relocate(Obesity, .after=Diabetes) %>%
  relocate(Hearing_impairment, .after=Obesity) %>%
  relocate(Smoking, .after=Hearing_impairment) %>%
  relocate(Depression, .after=Smoking) %>%
  select(-c(Systolic_blood_pressure, Diastolic_blood_pressure, total_cholesterol_value, 
            LDL_value, Height, Weight, BMI, hba1c_value, treatment_diabetes, treatment_hypertension,
            past_depression, smoking_present, epoch_score_currently, pittsburgh_total_score,
            Sleep_disturbances_holder, Cognitive_activity_holder, Social_activity_total,
            Social_activity_holder, social_life_frequency_activities:social_life_frequency_phone_calls,
            diagnosed_impairment, treatment_hyperlipidemia, gds_score))


## Recoding exercise columns
exercise_libra <- lifestyle_raw %>%
  select(CONP_ID, exer_curr_act1_intensity:exer_curr_act5_hours) %>%
  mutate(
    moderate_minutes_week = 0,
    vigorous_minutes_week = 0
  ) %>%
  mutate(
    across(contains("_intensity"), ~ as.numeric(as.character(.))),
    across(contains("_days"), ~ as.numeric(as.character(.))),
    across(contains("_hours"), ~ as.numeric(as.character(.)))
  ) %>%
  rowwise() %>%
  mutate(
    has_any_exercise_data = any(!is.na(c_across(contains("exer_curr_act"))))
  ) %>%
  ungroup()

# Loop through activities to calculate moderate and vigorous minutes
for(i in 1:5) {  
  intensity_col <- paste0("exer_curr_act", i, "_intensity")
  days_col <- paste0("exer_curr_act", i, "_days")
  hours_col <- paste0("exer_curr_act", i, "_hours")
  
  exercise_libra <- exercise_libra %>%
    mutate(
      # Calculate minutes per week for this activity
      activity_minutes = case_when(
        !is.na(.data[[intensity_col]]) & 
          !is.na(.data[[days_col]]) & 
          !is.na(.data[[hours_col]]) ~ .data[[days_col]] * .data[[hours_col]] * 60,
        TRUE ~ 0
      ),
      
      # Add to moderate if intensity == 2
      moderate_minutes_week = moderate_minutes_week + 
        ifelse(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 2, 
               activity_minutes, 0),
      
      # Add to vigorous if intensity == 3
      vigorous_minutes_week = vigorous_minutes_week + 
        ifelse(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 3, 
               activity_minutes, 0)
    ) %>%
    select(-activity_minutes)
}

# LIBRA2 exercise criteria
exercise_libra <- exercise_libra %>%
  mutate(
    # Check if meets any guideline:
    # 1) >= 120 min/week moderate (2 hr per week)
    # 2) >= 60 min/week vigorous (1 hr per week)  
    meets_moderate_guideline = moderate_minutes_week >= 150,
    meets_vigorous_guideline = vigorous_minutes_week >= 60,
    
    # LIBRA physical activity (1 = risk factor present, 0 = absent)
    Physical_inactivity = case_when(
      # If any guideline is met, no risk (0)
      meets_moderate_guideline | meets_vigorous_guideline ~ 0,
      # If no guideline is met, risk present (1)
      TRUE ~ 1
    )
  ) %>%
  select(CONP_ID, Physical_inactivity)

# Merge
dat_libra <- dat_libra %>%
  left_join(exercise_libra, by="CONP_ID")

# Save
saveRDS(dat_libra, 
        file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_libra_dat.rds"))
