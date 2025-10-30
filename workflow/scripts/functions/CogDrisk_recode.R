
# CogDrisk ----------------------------------------------------------------
## Recode dataframe according to CogDrisk values


# PREVENT-AD COGDRISK Data Preparation and Recoding Functions
# Modular approach for flexible imputation testing

# ============================================================================
# 1. CLINICAL DATA RECODING FUNCTION
# ============================================================================

recode_clinical_cogdrisk <- function(clinical_raw) {
  clinical_recoded <- clinical_raw %>%
    mutate(
      # Sex: 0 = Female, 1 = Male
      Sex = if_else(Sex == "Female", 0, 1),
      
      # Age: Convert from months to years
      Age = Candidate_Age / 12,
      
      # Education level (categorical)
      Education_level = case_when(
        Education_years > 11 ~ "High",
        Education_years > 8 & Education_years <= 11 ~ "Middle",
        Education_years <= 8 ~ "Low",
        TRUE ~ NA_character_
      ),
      
      # Hypertension: Based on blood pressure thresholds
      Hypertension = if_else(
        Systolic_blood_pressure > 129 | Diastolic_blood_pressure >= 80,
        1, 0
      ),
      
      # BMI calculation
      BMI = Weight / (Height / 100)^2,
      
      # BMI category
      BMI_category = case_when(
        BMI < 18.5 ~ "Underweight",
        BMI >= 18.5 & BMI < 25 ~ "Normal",
        BMI >= 25 & BMI < 30 ~ "Overweight",
        BMI >= 30 ~ "Obese",
        TRUE ~ NA_character_
      ),
      
      # High cholesterol: >6.5 mmol/L
      High_cholesterol = if_else(total_cholesterol_value > 6.5, 1, 0),
      
      # Depression (from past history)
      Depression = past_depression,
      
      # Atrial fibrillation (from past history)
      Atrial_fibrillation = past_atrial_fibrillation,
      
      # Diabetes: Based on treatment or HbA1c by age
      Diabetes = case_when(
        treatment_diabetes > 0 ~ 1,
        Age < 60 & hba1c_value <= 0.061 ~ 0,
        (Age >= 60 & Age < 70) & hba1c_value <= 0.075 ~ 0,
        Age >= 70 & hba1c_value <= 0.070 ~ 0,
        TRUE ~ 1
      ),
      
      # TBI: Traumatic brain injury
      TBI = if_else(
        head_injury_hospitalized == 1 | head_injury_severe == 1,
        1, 0
      )
    ) %>%
    # Reorganize columns
    relocate(Age, .after = Sex) %>%
    relocate(Education_level, .after = Education_years) %>%
    relocate(BMI, .after = Weight) %>%
    relocate(BMI_category, .after = BMI) %>%
    relocate(Hypertension, .after = BMI_category) %>%
    relocate(High_cholesterol, .after = total_cholesterol_value) %>%
    relocate(Depression, .after = past_depression) %>%
    relocate(Atrial_fibrillation, .after = past_atrial_fibrillation) %>%
    relocate(Diabetes, .after = treatment_diabetes) %>%
    # Remove raw/intermediate columns
    select(-c(
      Candidate_Age, Systolic_blood_pressure, Diastolic_blood_pressure,
      Height, Weight, total_cholesterol_value, past_depression,
      past_atrial_fibrillation, treatment_diabetes, hba1c_value,
      head_injury_hospitalized, head_injury_severe
    ))
  
  return(clinical_recoded)
}


# ============================================================================
# 2. LIFESTYLE DATA RECODING FUNCTION
# ============================================================================

recode_lifestyle_cogdrisk <- function(lifestyle_raw) {
  
  # Smoking recoding
  lifestyle_recoded <- lifestyle_raw %>%
    mutate(
      # Smoking: 0 = non-smoker, 1 = former, 2 = current
      Smoking = case_when(
        smoking_present == 4 ~ 2,  # Current
        smoking_present == 3 ~ 2,  # Current
        smoking_present == 2 ~ 1,  # Former
        smoking_present == 1 ~ 1,  # Former
        smoking_present == 0 ~ 0,  # Non-smoker
        TRUE ~ NA_real_
      ),
      
      # Cognitive engagement: Tertiles
      Cognitive_engagement = ntile(epoch_score_currently, 3),
      
      # Social engagement: Calculate mean of social activities
      Social_engagement_holder = rowMeans(
        select(., social_life_frequency_activities:social_life_frequency_phone_calls),
        na.rm = TRUE
      ),
      
      # Social engagement: 0 = not lonely, 1 = lonely
      Social_engagement = if_else(
        Social_engagement_holder > mean(Social_engagement_holder, na.rm = TRUE),
        0, 1
      )
    ) %>%
    select(-Social_engagement_holder)
  
  # Initialize exercise columns
  lifestyle_recoded <- lifestyle_recoded %>%
    mutate(
      light_minutes_week = 0,
      moderate_minutes_week = 0,
      heavy_minutes_week = 0
    ) %>%
    mutate(
      across(contains("_intensity"), ~ as.numeric(as.character(.))),
      across(contains("_days"), ~ as.numeric(as.character(.))),
      across(contains("_hours"), ~ as.numeric(as.character(.)))
    )
  
  # Calculate exercise minutes for each activity (1-5)
  for (i in 1:5) {
    intensity_col <- paste0("exer_curr_act", i, "_intensity")
    days_col <- paste0("exer_curr_act", i, "_days")
    hours_col <- paste0("exer_curr_act", i, "_hours")
    
    lifestyle_recoded <- lifestyle_recoded %>%
      mutate(
        # Calculate minutes only if all three values are present
        activity_minutes = case_when(
          !is.na(.data[[intensity_col]]) &
            !is.na(.data[[days_col]]) &
            !is.na(.data[[hours_col]]) ~ .data[[days_col]] * .data[[hours_col]] * 60,
          TRUE ~ 0
        ),
        
        # Add to appropriate intensity category
        light_minutes_week = light_minutes_week +
          if_else(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 1,
                  activity_minutes, 0),
        
        moderate_minutes_week = moderate_minutes_week +
          if_else(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 2,
                  activity_minutes, 0),
        
        heavy_minutes_week = heavy_minutes_week +
          if_else(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 3,
                  activity_minutes, 0)
      ) %>%
      select(-activity_minutes)
  }
  
  # Final processing
  lifestyle_recoded <- lifestyle_recoded %>%
    relocate(Smoking, .after = smoking_present) %>%
    relocate(Cognitive_engagement, .after = epoch_score_currently) %>%
    relocate(light_minutes_week, .after = Cognitive_engagement) %>%
    relocate(moderate_minutes_week, .after = light_minutes_week) %>%
    relocate(heavy_minutes_week, .after = moderate_minutes_week) %>%
    # Physical inactivity: 0 = active (>150 min/week), 1 = inactive
    mutate(
      Physical_inactivity = if_else(
        (moderate_minutes_week + heavy_minutes_week) > 150,
        0, 1
      )
    ) %>%
    # Keep only relevant columns
    select(CONP_ID, Smoking, Cognitive_engagement, Social_engagement, Physical_inactivity)
  
  return(lifestyle_recoded)
}


# ============================================================================
# 3. MASTER RECODING FUNCTION
# ============================================================================

prepare_preventad_cogdrisk <- function(clinical_raw, lifestyle_raw, verbose = TRUE) {
  
  if (verbose) {
    message("=== Recoding Clinical Data ===")
  }
  clinical_recoded <- recode_clinical_cogdrisk(clinical_raw)
  
  if (verbose) {
    message("=== Recoding Lifestyle Data ===")
  }
  lifestyle_recoded <- recode_lifestyle_cogdrisk(lifestyle_raw)
  
  if (verbose) {
    message("=== Merging Datasets ===")
  }
  
  # Merge clinical and lifestyle data
  data_cogdrisk <- clinical_recoded %>%
    left_join(lifestyle_recoded, by = "CONP_ID")
  
  if (verbose) {
    message("\n=== Data Preparation Complete ===")
    message("Total rows: ", nrow(data_cogdrisk))
    message("Total columns: ", ncol(data_cogdrisk))
    
    # Report missingness
    message("\n=== Missingness Summary ===")
    missing_summary <- colSums(is.na(data_cogdrisk))
    missing_summary <- missing_summary[missing_summary > 0]
    
    if (length(missing_summary) > 0) {
      for (var in names(missing_summary)) {
        pct <- round(missing_summary[var] / nrow(data_cogdrisk) * 100, 1)
        message(sprintf("  %s: %d (%.1f%%)", var, missing_summary[var], pct))
      }
    } else {
      message("  No missing values detected!")
    }
  }
  
  return(data_cogdrisk)
}
