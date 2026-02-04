
# SETUP -------------------------------------------------------------------

# Load files
source('workflow/scripts/config.R')
source(CRS_FN)

# Load imputed CRS factors (contains all clinical, lifestyle, and medication data)
crs_factors <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_crsfactors_imputed.rds"))


# RECODING ----------------------------------------------------------------
dat_libra <- crs_factors %>%
  mutate(
    # Age = age,
    # Sleep disturbances:
    # 1 == highest tertile
    Sleep_disturbances_holder = ntile(pittsburgh_total_score, 3),
    Sleep_disturbances = if_else(Sleep_disturbances_holder == 3, 1, 0),
    
    # Hypertension risk factor present == 1
    Hypertension = case_when(
      Systolic_blood_pressure >= 140 ~ 1,
      Diastolic_blood_pressure >= 90 ~ 1,
      treatment_hypertension == 1 | treatment_hypertension == 2 ~ 1,
      medusage_antihypertensive == 1 ~ 1,
      is.na(Systolic_blood_pressure) & is.na(Diastolic_blood_pressure) & is.na(treatment_hypertension) & is.na(medusage_antihypertensive) ~ NA,
      TRUE ~ 0
    ),
    
    # Hypercholesterolemia
    # 1 == TC >= 6.2 or currently treating hyperlipidemia
    Hypercholesterolemia = case_when(
      total_cholesterol_value >= 6.2 ~ 1,
      treatment_hyperlipidemia == 1 | treatment_hyperlipidemia == 2 ~ 1,
      medusage_antihyperlipidemic == 1 ~ 1,
      is.na(total_cholesterol_value) & is.na(treatment_hyperlipidemia) & is.na(medusage_antihyperlipidemic) ~ NA,
      TRUE ~ 0
    ),
    
    BMI = Weight/(Height/100)^2,
    Obesity = if_else(BMI >= 30, 1, 0),
    
    # Hearing impairment:
    # 1 == yes 
    Hearing_impairment = if_else(diagnosed_impairment != 0, 1, 0),
    
    Diabetes = case_when(
      treatment_diabetes > 0 ~ 1,
      medusage_antidiabetic == 1 ~ 1,
      Age < 60 & treatment_diabetes == 0 & hba1c_value <= 0.061 ~ 0,
      (Age >= 60 & Age < 70) & hba1c_value <= 0.075 & treatment_diabetes == 0 ~ 0,
      Age >= 70 & treatment_diabetes == 0 & hba1c_value <= 0.07 ~ 0,
      is.na(treatment_diabetes) & is.na(medusage_antidiabetic) & is.na(hba1c_value) ~ NA,
      TRUE ~ 1),
    
    # Depression
    # 1 == yes
    Depression = case_when(past_depression == 1 ~ 1,
                           medusage_antidepressants == 1 ~ 1,
                           gds_score >= 5 ~ 1,
                           past_depression == 0 & medusage_antidepressants == 0 ~ 0,
                           TRUE ~ NA),
    
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
         Social_participation = if_else(Social_activity_holder == 1, 1, 0)) %>%
  
  relocate(Sleep_disturbances, .after=CONP_ID) %>%
  relocate(Hypertension, .after=Diastolic_blood_pressure) %>%
  relocate(Hypercholesterolemia, .after=LDL_value) %>%
  relocate(BMI, .after=Hypercholesterolemia) %>%
  relocate(Diabetes, .after=BMI) %>%
  relocate(Obesity, .after=Diabetes) %>%
  relocate(Hearing_impairment, .after=Obesity) %>%
  relocate(Smoking, .after=Hearing_impairment) %>%
  relocate(Depression, .after=Smoking) %>%
  select(CONP_ID, Sleep_disturbances, Hypertension, Hypercholesterolemia, Obesity, 
         Hearing_impairment, Diabetes, Depression, Smoking, Cognitive_activity, Social_participation)


## Recoding exercise columns
exercise_libra <- crs_factors %>%
  select(CONP_ID, exer_curr_act1_intensity:exer_curr_act5_hours) %>%
  mutate(
    moderate_minutes_week = 0,
    vigorous_minutes_week = 0
  ) %>%
  # Ensure all exercise variables are numeric (should already be from crs_factors_imputed)
  mutate(
    across(contains("_intensity"), ~ as.numeric(.)),
    across(contains("_weeks"), ~ as.numeric(.)),
    across(contains("_hours"), ~ as.numeric(.))
  ) %>%
  rowwise() %>%
  mutate(
    has_any_exercise_data = any(!is.na(c_across(contains("exer_curr_act"))))
  ) %>%
  ungroup()

# Loop through activities to calculate moderate and vigorous minutes
for(i in 1:5) {
  intensity_col <- paste0("exer_curr_act", i, "_intensity")
  weeks_col <- paste0("exer_curr_act", i, "_weeks")
  hours_col <- paste0("exer_curr_act", i, "_hours")
  
  exercise_libra <- exercise_libra %>%
    mutate(
      # Calculate minutes per week for this activity
      activity_minutes = case_when(
        !is.na(.data[[intensity_col]]) &
          !is.na(.data[[weeks_col]]) &
          !is.na(.data[[hours_col]]) ~ .data[[weeks_col]] * .data[[hours_col]] * 60,
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
    Physical_activity = case_when(
      # If any guideline is met, no risk (0)
      meets_moderate_guideline | meets_vigorous_guideline ~ 0,
      # If no guideline is met, risk present (1)
      TRUE ~ 1
    )
  ) %>%
  select(CONP_ID, Physical_activity)

# Merge
dat_libra <- dat_libra %>%
  left_join(exercise_libra, by="CONP_ID")

# Scoring
dat_libra_scored <- dat_libra %>%
  mutate(score_libra2 = calculate_libra2(.)) %>%
  relocate(score_libra2, .after=CONP_ID) 

# Save
saveRDS(dat_libra, 
        file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_libra2.rds"))
saveRDS(dat_libra_scored, 
        file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_libra2_scored.rds"))
