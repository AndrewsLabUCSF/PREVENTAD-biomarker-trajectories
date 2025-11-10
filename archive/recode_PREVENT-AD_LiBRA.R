
# SETUP -------------------------------------------------------------------

# Load files
source('workflow/scripts/config.R')

clinical_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_clinical_raw.rds"))
lifestyle_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_lifestyle_raw.rds"))


# RECODING ----------------------------------------------------------------
dat_libra <- clinical_raw %>%
  select(CONP_ID, Systolic_blood_pressure, Diastolic_blood_pressure, total_cholesterol_value, 
         LDL_value, Height, Weight, hba1c_value, treatment_diabetes, past_depression) %>%
  left_join((lifestyle_raw %>%
               select(CONP_ID, smoking_present, epoch_score_currently)),
            by="CONP_ID") %>%
  mutate(
    Hypertension = if_else(
      Systolic_blood_pressure >= 140 | Diastolic_blood_pressure >= 90, 1, 0),
    Hypercholesterolemia = if_else(
      total_cholesterol_value >= 5.2 & LDL_value >= 4.2, 1, 0),
    BMI = Weight/(Height/100)^2,
    Obesity = if_else(BMI >= 30, 1, 0),
    Diabetes = case_when(treatment_diabetes > 0 ~ 1,
                         hba1c_value >= 0.065 ~ 1,
                         TRUE ~ 0),
    Depression = past_depression %>% factor(),
    Smoking = if_else(smoking_present >= 3, 1, 0),
    Cognitive_activity_holder = ntile(epoch_score_currently, 3),
    
    # Cognitive_activity:
    # 1 == highest tertile
    Cognitive_activity = if_else(Cognitive_activity_holder == 3, 1, 0)
  ) %>%
  relocate(Hypertension, .after=Diastolic_blood_pressure) %>%
  relocate(Hypercholesterolemia, .after=LDL_value) %>%
  relocate(BMI, .after=Weight) %>%
  relocate(Obesity, .after=BMI) %>%
  relocate(Diabetes, .after=Obesity) %>%
  relocate(Depression, .after=past_depression) %>%
  relocate(Smoking, .after=smoking_present) %>%
  select(-c(Systolic_blood_pressure, Diastolic_blood_pressure, total_cholesterol_value, 
            LDL_value, Height, Weight, BMI, hba1c_value, treatment_diabetes,
            past_depression, smoking_present, epoch_score_currently, 
            Cognitive_activity_holder))


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

# Apply LIBRA criteria
exercise_libra <- exercise_libra %>%
  mutate(
    # Check if meets any guideline:
    # 1) ≥150 min/week moderate (30 min × 5 days)
    # 2) ≥60 min/week vigorous (20 min × 3 days)  
    # 3) Combination: moderate + (vigorous × 2) ≥ 150
    #    (vigorous counts double toward moderate equivalent)
    meets_moderate_guideline = moderate_minutes_week >= 150,
    meets_vigorous_guideline = vigorous_minutes_week >= 60,
    meets_combined_guideline = (moderate_minutes_week + (vigorous_minutes_week * 2)) >= 150,
    
    # LIBRA physical inactivity (1 = risk factor present, 0 = absent)
    Physical_inactivity = case_when(
      # If any guideline is met, no risk (0)
      meets_moderate_guideline | meets_vigorous_guideline | meets_combined_guideline ~ 0,
      # If no guideline is met, risk present (1)
      TRUE ~ 1
    )
  ) %>%
  select(CONP_ID, Physical_inactivity)

exercise_libra <- exercise_libra %>%
  mutate(
    light_minutes_week = if_else(has_any_exercise_data, light_minutes_week, NA_real_),
    moderate_minutes_week = if_else(has_any_exercise_data, moderate_minutes_week, NA_real_),
    heavy_minutes_week = if_else(has_any_exercise_data, heavy_minutes_week, NA_real_)
  ) %>%
  select(-has_any_exercise_data)

# Merge
dat_libra <- dat_libra %>%
  left_join(exercise_libra, by="CONP_ID")

# Save
saveRDS(dat_libra, 
        file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_libra_dat.rds"))
