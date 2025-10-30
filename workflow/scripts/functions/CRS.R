
# CogDrisk ----------------------------------------------------------------
## Calculate CogDrisk score
##  Inputs:

calculate_cogdrisk <- function(df, 
                               age_var="Age",
                               sex_var="Sex",
                               educ_var="Education_level",
                               obesity_var="BMI_category",
                               chol_var="High_cholesterol",
                               diabetes_var="Diabetes",
                               stroke_var="Stroke",
                               tbi_var="TBI",
                               hypertension_var="Hypertension",
                               atrialfib_var="Atrial_fibrillation",
                               insomnia_var="Insomnia",
                               depression_var="Depression",
                               physact_var="Physical_inactivity",
                               cogeng_var="Cognitive_engagement",
                               soceng_var="Social_engagement",
                               diet_var="Diet",
                               smoking_var="Smoking",
                               verbose=TRUE
                               ) {
  
  required_vars <- c(age_var, sex_var, educ_var, obesity_var, chol_var,
                     diabetes_var, stroke_var, tbi_var, hypertension_var,
                     atrialfib_var, insomnia_var, depression_var, physact_var, 
                     cogeng_var, soceng_var, diet_var, smoking_var)
  
  var_names <- c("Age", "Sex", "Education", "Obesity",
                 "High cholesterol", "Diabetes", "Stroke", "TBI", "Hypertension",
                 "Atrial fibrillation", "Insomnia", "Depression", "Physical inactivity", 
                 "Cognitive engagement", "Social engagement", "Diet", "Smoking")
  
  available_vars <- required_vars[required_vars %in% names(df)]
  missing_vars <- setdiff(required_vars, names(df))
  
  # Warning message for missing variables
  if (length(missing_vars) > 0 && verbose) {
    missing_names <- var_names[match(missing_vars, required_vars)]
    message("Calculating partial COGDRISK score.")
    message("Missing variables: ", paste(missing_names, collapse = ", "))
    message("Score will be calculated using ", length(available_vars), 
            " out of ", length(required_vars), " variables.\n")
  }
  
  # Initialize points
  points <- rep(0, nrow(df))
  components_used <- c()
  components_missing <- c()
  
  # Calculate score
  ## Age and Sex
  if (age_var %in% names(df) && sex_var %in% names(df)) {
    age <- df[[age_var]]
    sex <- df[[sex_var]]
    
    agesex_points <- case_when(
      # Female (0)
      sex == 0 & age < 65 ~ 0,
      sex == 0 & age >= 65 & age < 69 ~ 5,
      sex == 0 & age >= 70 & age < 74 ~ 7, 
      sex == 0 & age >= 75 & age < 79 ~ 13,
      sex == 0 & age >= 80 & age < 84 ~ 16,
      sex == 0 & age >= 85 & age < 89 ~ 19,
      sex == 0 & age >= 90 ~ 23,
      
      # Male (1)
      sex == 1 & age < 65 ~ 0,
      sex == 1 & age >= 65 & age < 69 ~ 5,
      sex == 1 & age >= 70 & age < 74 ~ 8, 
      sex == 1 & age >= 75 & age < 79 ~ 12,
      sex == 1 & age >= 80 & age < 84 ~ 17,
      sex == 1 & age >= 85 & age < 89 ~ 19,
      sex == 1 & age >= 90 ~ 23,
      
      TRUE ~ 0
    )
    points <- points + agesex_points
    components_used <- c(components_used, "age", "sex")
  } else {
    components_missing <- c(components_missing, "age", "sex")
  }
  
  ## Education
  if (educ_var %in% names(df)) {
    educ <- df[[educ_var]]
    educ_points <- case_when(
      educ == "High" ~ 0,
      educ == "Middle" ~ 2,
      educ == "Low" ~ 4,
      TRUE ~ 0
    )
    points <- points + educ_points
    components_used <- c(components_used, "education")
  } else {
    components_missing <- c(components_missing, "education")
  }
  
  ## Obesity
  if (obesity_var %in% names(df)) {
    obesity <- df[[obesity_var]]
    obesity_points <- case_when(
      obesity == "Normal" ~ 0,
      obesity == "Overweight" ~ 1,
      obesity == "Underweight" ~ 3,
      obesity == "Obese" ~ 2,
      TRUE ~ 0
    )
    points <- points + obesity_points
    components_used <- c(components_used, "obesity")
  } else {
    components_missing <- c(components_missing, "obesity")
  }
  
  ## High cholesterol
  if (chol_var %in% names(df)) {
    chol <- df[[chol_var]]
    chol_points <- if_else(chol == 0, 0, 1, missing=0)
    points <- points + chol_points
    components_used <- c(components_used, "")
  }
  
  ## Diabetes
  if (diabetes_var %in% names(df)) {
    diabetes <- df[[diabetes_var]]
    diabetes_points <- case_when(
      sex == 0 & diabetes == 1 ~ 3,
      sex == 1 & diabetes == 1 ~ 2,
      diabetes == 0 ~ 0,
      TRUE ~ 0
    )
    points <- points + diabetes_points
    components_used <- c(components_used, "diabetes")
  } else {
    components_missing <- c(components_missing, "diabetes")
  }
  
  ## Stroke
  
  ## TBI
  if (tbi_var %in% names(df)) {
    tbi <- df[[tbi_var]]
    tbi_points <- if_else(tbi == 0, 0, 1, missing=0)
    points <- points + tbi_points
    components_used <- c(components_used, "TBI")
  } else {
    components_missing <- c(components_missing, "TBI")
  }
  
  ## Hypertension
  if (hypertension_var %in% names(df)) {
    hypertension <- df[[hypertension_var]]
    hypertension_points <- if_else(hypertension == 0, 0, 1, missing=0)
    points <- points + hypertension_points
    components_used <- c(components_used, "hypertension")
  } else {
    components_missing <- c(components_missing, "hypertension")
  }
  
  ## Atrial fibrillation
  if (atrialfib_var %in% names(df)) {
    atrialfib <- df[[atrialfib_var]]
    atrialfib_points <- if_else(atrialfib == 0, 0, 1, missing=0)
    points <- points + atrialfib_points
    components_used <- c(components_used, "atrial fibrillation")
  } else {
    components_missing <- c(components_missing, "atrial fibrillation")
  }
  
  ## Insomnia
  
  ## Depression
  if (depression_var %in% names(df)) {
    depression <- df[[depression_var]]
    depression_points <- if_else(depression == 0, 0, 3, missing=0)
    points <- points + depression_points
    components_used <- c(components_used, "depression")
  } else {
    components_missing <- c(components_missing, "depression")
  }
  
  ## Physical inactivity
  if (physact_var %in% names(df)) {
    physact <- df[[physact_var]]
    physact_points <- if_else(physact == 0, 0, -3, missing=0)
    points <- points + physact_points
    components_used <- c(components_used, "physical activity")
  } else {
    components_missing <- c(components_missing, "physical activity")
  }
  
  ## Cognitive engagement
  if (cogeng_var %in% names(df)) {
    cogeng <- df[[cogeng_var]]
    cogeng_points <- case_when(
      cogeng == 1 ~ 0,
      cogeng == 2 ~ -5,
      cogeng == 3 ~ -4,
      TRUE ~ 0
    )
    points <- points + cogeng_points
    components_used <- c(components_used, "cognitive engagement")
  } else {
    components_missing <- c(components_missing, "cognitive engagement")
  }
  
  ## Social engagement
  if (soceng_var %in% names(df)) {
    soceng <- df[[soceng_var]]
    soceng_points <- if_else(soceng == 0, 0, 2, missing=0)
    points <- points + soceng_points
    components_used <- c(components_used, "social engagement")
  } else {
    components_missing <- c(components_missing, "social engagement")
  }
  
  ## Diet
  
  ## Smoking
  if (smoking_var %in% names(df)) {
    smoking <- df[[smoking_var]]
    smoking_points <- case_when(
      smoking == 2 ~ 1,
      smoking == 1 ~ 0,
      smoking == 0 ~ 0,
      TRUE ~ 0
    )
    points <- points + smoking_points
    components_used <- c(components_used, "smoking")
  } else {
    components_missing <- c(components_missing, "smoking")
  }
  
  return(points)

}

# test <- data_cogdrisk %>% 
#   mutate(score = calculate_cogdrisk(data_cogdrisk)) %>%
#   relocate(score, .after=CONP_ID)


# LIBRA -------------------------------------------------------------------
calculate_libra <- function(df, 
                            hypertension_var="Hypertension",
                            hyperchol_var="Hypercholesterolemia",
                            obesity_var="Obesity",       
                            physact_var="Physical_inactivity",
                            diabetes_var="Diabetes",
                            depression_var="Depression",
                            smoking_var="Smoking",
                            alcohol_var="Alcohol",
                            cogact_var="Cognitive_activity",
                            diet_var="Diet",
                            chd_var="Chd",
                            renaldys_var="Renal_dysfunction",
                            verbose=TRUE
                            ) {
  
  required_vars <- c(hypertension_var, hyperchol_var, obesity_var, physact_var, 
                     diabetes_var, depression_var, smoking_var, alcohol_var,
                     cogact_var, diet_var, chd_var, renaldys_var)
  
  var_names <- c("Hypertension", "Hypercholesterolemia", "Obesity", "Physical inactivity", 
                 "Diabetes", "Depression", "Smoking", "Alcohol", "Cognitive activity",
                 "Diet", "Chd", "Renal dysfunction")
  
  available_vars <- required_vars[required_vars %in% names(df)]
  missing_vars <- setdiff(required_vars, names(df))
  
  # Warning message for missing variables
  if (length(missing_vars) > 0 && verbose) {
    missing_names <- var_names[match(missing_vars, required_vars)]
    message("Calculating partial LIBRA score.")
    message("Missing variables: ", paste(missing_names, collapse = ", "))
    message("Score will be calculated using ", length(available_vars), 
            " out of ", length(required_vars), " variables.\n")
  }
  
  # Initialize points
  points <- rep(0, nrow(df))
  components_used <- c()
  components_missing <- c()
  
  # Calculate score
  ## Hypertension
  if (hypertension_var %in% names(df)) {
    hypertension <- df[[hypertension_var]]
    hypertension_points <- if_else(hypertension == 0, 0, 1.6, missing=0)
    points <- points + hypertension_points
    components_used <- c(components_used, "hypertension")
  } else {
    components_missing <- c(components_missing, "hypertension")
  }
  
  ## Hypercholesterolemia
  if (hyperchol_var %in% names(df)) {
    hyperchol <- df[[hyperchol_var]]
    hyperchol_points <- if_else(hyperchol == 0, 0, 1.4, missing=0)
    points <- points + hyperchol_points
    components_used <- c(components_used, "hypercholesterolemia")
  }
  
  ## Obesity
  if (obesity_var %in% names(df)) {
    obesity <- df[[obesity_var]]
    obesity_points <- if_else(obesity == 0, 0, 1.3, missing=0)
    points <- points + obesity_points
    components_used <- c(components_used, "obesity")
  } else {
    components_missing <- c(components_missing, "obesity")
  }
  
  ## Physical inactivity
  if (physact_var %in% names(df)) {
    physact <- df[[physact_var]]
    physact_points <- if_else(physact == 0, 0, 1.1, missing=0)
    points <- points + physact_points
    components_used <- c(components_used, "physical activity")
  } else {
    components_missing <- c(components_missing, "physical activity")
  }
  
  ## Diabetes
  if (diabetes_var %in% names(df)) {
    diabetes <- df[[diabetes_var]]
    diabetes_points <- if_else(diabetes == 0, 0, 1.3, missing=0)
    points <- points + diabetes_points
    components_used <- c(components_used, "diabetes")
  } else {
    components_missing <- c(components_missing, "diabetes")
  }
  
  ## Depression
  if (depression_var %in% names(df)) {
    depression <- df[[depression_var]]
    depression_points <- if_else(depression == 0, 0, 2.1, missing=0)
    points <- points + depression_points
    components_used <- c(components_used, "depression")
  } else {
    components_missing <- c(components_missing, "depression")
  }
  
  ## Smoking
  if (smoking_var %in% names(df)) {
    smoking <- df[[smoking_var]]
    smoking_points <- if_else(smoking == 0, 0, 1.5, missing=0)
    points <- points + smoking_points
    components_used <- c(components_used, "smoking")
  } else {
    components_missing <- c(components_missing, "smoking")
  }
  
  ## Alcohol 
  
  ## Cognitive activity
  if (cogact_var %in% names(df)) {
    cogact <- df[[cogact_var]]
    cogact_points <- if_else(cogact == 0, 0, -3.2, missing=0)
    points <- points + cogact_points
    components_used <- c(components_used, "cognitive activity")
  } else {
    components_missing <- c(components_missing, "cognitive activity")
  }
  
  ## Diet
  
  ## CHD
  
  ## Renal dysfunction
  
  return(points)
  
}

# test <- dat_libra %>%
#   mutate(score = calculate_libra(dat_libra)) %>%
#   relocate(score, .after=CONP_ID)
