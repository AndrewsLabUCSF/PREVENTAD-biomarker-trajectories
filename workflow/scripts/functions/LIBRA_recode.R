# PREVENT-AD LIBRA Data Preparation and Recoding Functions
# Modular approach for flexible imputation testing

# ============================================================================
# 1. CLINICAL/LIFESTYLE DATA RECODING FUNCTION
# ============================================================================

recode_libra_clinical <- function(clinical_raw, lifestyle_raw) {
  
  dat_libra <- clinical_raw %>%
    select(CONP_ID, Systolic_blood_pressure, Diastolic_blood_pressure, 
           total_cholesterol_value, LDL_value, Height, Weight, hba1c_value, 
           treatment_diabetes, past_depression) %>%
    left_join(
      lifestyle_raw %>%
        select(CONP_ID, smoking_present, epoch_score_currently),
      by = "CONP_ID"
    ) %>%
    mutate(
      # Hypertension: SBP >= 140 or DBP >= 90
      Hypertension = if_else(
        Systolic_blood_pressure >= 140 | Diastolic_blood_pressure >= 90, 
        1, 0
      ),
      
      # Hypercholesterolemia: Total chol >= 5.2 AND LDL >= 4.2
      Hypercholesterolemia = if_else(
        total_cholesterol_value >= 5.2 & LDL_value >= 4.2, 
        1, 0
      ),
      
      # BMI calculation
      BMI = Weight / (Height / 100)^2,
      
      # Obesity: BMI >= 30
      Obesity = if_else(BMI >= 30, 1, 0),
      
      # Diabetes: Treatment or HbA1c >= 6.5%
      Diabetes = case_when(
        treatment_diabetes > 0 ~ 1,
        hba1c_value >= 0.065 ~ 1,
        TRUE ~ 0
      ),
      
      # Depression
      Depression = past_depression,
      
      # Smoking: Current smoker (3-4)
      Smoking = if_else(smoking_present >= 3, 1, 0),
      
      # Cognitive activity tertiles
      Cognitive_activity_holder = ntile(epoch_score_currently, 3),
      
      # Cognitive activity: 1 = highest tertile (protective)
      Cognitive_activity = if_else(Cognitive_activity_holder == 3, 1, 0)
    ) %>%
    # Reorganize columns
    relocate(Hypertension, .after = Diastolic_blood_pressure) %>%
    relocate(Hypercholesterolemia, .after = LDL_value) %>%
    relocate(BMI, .after = Weight) %>%
    relocate(Obesity, .after = BMI) %>%
    relocate(Diabetes, .after = Obesity) %>%
    relocate(Depression, .after = past_depression) %>%
    relocate(Smoking, .after = smoking_present) %>%
    # Remove raw/intermediate columns
    select(-c(
      Systolic_blood_pressure, Diastolic_blood_pressure, 
      total_cholesterol_value, LDL_value, Height, Weight, BMI, 
      hba1c_value, treatment_diabetes, past_depression, 
      smoking_present, epoch_score_currently, Cognitive_activity_holder
    ))
  
  return(dat_libra)
}


# ============================================================================
# 2. EXERCISE DATA RECODING FUNCTION
# ============================================================================

recode_libra_exercise <- function(lifestyle_raw) {
  
  # Select and initialize exercise columns
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
    )
  
  # Loop through activities to calculate moderate and vigorous minutes
  for (i in 1:5) {
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
          if_else(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 2,
                  activity_minutes, 0),
        
        # Add to vigorous if intensity == 3
        vigorous_minutes_week = vigorous_minutes_week +
          if_else(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 3,
                  activity_minutes, 0)
      ) %>%
      select(-activity_minutes)
  }
  
  # Apply LIBRA physical activity criteria
  exercise_libra <- exercise_libra %>%
    mutate(
      # Check if meets any guideline:
      # 1) >= 150 min/week moderate (30 min × 5 days)
      # 2) >= 60 min/week vigorous (20 min × 3 days)
      # 3) Combination: moderate + (vigorous × 2) >= 150
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
  
  return(exercise_libra)
}


# ============================================================================
# 3. MASTER RECODING FUNCTION
# ============================================================================

prepare_preventad_libra <- function(clinical_raw, lifestyle_raw, verbose = TRUE) {
  
  if (verbose) {
    message("=== Recoding LIBRA Clinical/Lifestyle Data ===")
  }
  dat_libra_clinical <- recode_libra_clinical(clinical_raw, lifestyle_raw)
  
  if (verbose) {
    message("=== Recoding LIBRA Exercise Data ===")
  }
  exercise_libra <- recode_libra_exercise(lifestyle_raw)
  
  if (verbose) {
    message("=== Merging Datasets ===")
  }
  
  # Merge clinical and exercise data
  dat_libra <- dat_libra_clinical %>%
    left_join(exercise_libra, by = "CONP_ID")
  
  if (verbose) {
    message("\n=== LIBRA Data Preparation Complete ===")
    message("Total rows: ", nrow(dat_libra))
    message("Total columns: ", ncol(dat_libra))
    
    # Report missingness
    message("\n=== Missingness Summary ===")
    missing_summary <- colSums(is.na(dat_libra))
    missing_summary <- missing_summary[missing_summary > 0]
    
    if (length(missing_summary) > 0) {
      for (var in names(missing_summary)) {
        pct <- round(missing_summary[var] / nrow(dat_libra) * 100, 1)
        message(sprintf("  %s: %d (%.1f%%)", var, missing_summary[var], pct))
      }
    } else {
      message("  No missing values detected!")
    }
    
    # Summary of risk factors
    message("\n=== LIBRA Risk Factor Prevalence ===")
    risk_factors <- c("Hypertension", "Hypercholesterolemia", "Obesity", 
                      "Diabetes", "Depression", "Smoking", "Physical_inactivity")
    
    for (rf in risk_factors) {
      if (rf %in% names(dat_libra)) {
        n_present <- sum(dat_libra[[rf]] == 1, na.rm = TRUE)
        n_total <- sum(!is.na(dat_libra[[rf]]))
        pct <- round(n_present / n_total * 100, 1)
        message(sprintf("  %s: %d/%d (%.1f%%)", rf, n_present, n_total, pct))
      }
    }
  }
  
  return(dat_libra)
}