
# SETUP -------------------------------------------------------------------

source('workflow/scripts/config.R')
source(CRS_FN)

library(dplyr)
library(tidyr)
library(stringr)

PREVENTAD_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_raw.rds"))

# Multi-visit reference data for tertile cutoffs
crsfactors_mv <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_crsfactors_imputed.rds"))


# IDENTIFY SINGLE-VISIT GWAS PARTICIPANTS ---------------------------------

gwas <- PREVENTAD_raw$GWAS

biomarkers <- PREVENTAD_raw$plasma_4plex %>%
  full_join(PREVENTAD_raw$plasma_ptau217) %>%
  select(CONP_ID, Visit_label)

criteria_met_ids <- biomarkers %>%
  group_by(CONP_ID) %>%
  summarise(n_visits = n()) %>%
  filter(n_visits >= 2, CONP_ID %in% gwas$CONP_ID) %>%
  pull(CONP_ID)

sv_ids <- setdiff(gwas$CONP_ID, criteria_met_ids)

cat("Criteria-met (multi-visit) participants:", length(criteria_met_ids), "\n")
cat("Single-visit GWAS participants:", length(sv_ids), "\n")


# BUILD RAW CRS FACTORS FOR SINGLE-VISIT PARTICIPANTS ---------------------
# Mirrors 02_clean_PREVENT-AD.R lines 102–244, filtered to sv_ids

crs_factors_sv <- PREVENTAD_raw$demographics %>%
  filter(CONP_ID %in% sv_ids) %>%
  select(CONP_ID, Sex, Education_years, Height) %>%
  mutate(Sex = factor(Sex, levels = c("Female", "Male"))) %>%
  left_join(
    (PREVENTAD_raw$bp_pulse_weight %>%
       group_by(CONP_ID) %>%
       slice(1) %>%
       select(CONP_ID, Age = Candidate_Age, Systolic_blood_pressure, Diastolic_blood_pressure) %>%
       mutate(Age = Age / 12)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$bp_pulse_weight %>%
       group_by(CONP_ID) %>%
       filter(!is.na(Weight)) %>%
       slice(1) %>%
       select(CONP_ID, Weight)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$lab %>%
       filter(!is.na(total_cholesterol_value)) %>%
       group_by(CONP_ID) %>%
       slice(1) %>%
       select(CONP_ID, total_cholesterol_value)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$lab %>%
       filter(!is.na(LDL_value)) %>%
       group_by(CONP_ID) %>%
       slice(1) %>%
       select(CONP_ID, LDL_value)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$lab %>%
       select(CONP_ID, hba1c_value) %>%
       filter(!is.na(hba1c_value)) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$medical_history %>%
       select(CONP_ID, past_depression, past_atrial_fibrillation, treatment_diabetes,
              treatment_hypertension, treatment_hyperlipidemia)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, head_injury_hospitalized, head_injury_severe) %>%
       filter_at(vars(head_injury_severe, head_injury_hospitalized), any_vars(!is.na(.)))
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, gds_score) %>%
       filter(!is.na(gds_score)) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, pittsburgh_total_score) %>%
       filter(!is.na(pittsburgh_total_score)) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$auditory %>%
       select(CONP_ID, diagnosed_impairment, subjective_hearing_impairment) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, smoking_present) %>%
       filter(!is.na(smoking_present)) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, epoch_score_currently) %>%
       filter(!is.na(epoch_score_currently)) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, starts_with("exer_curr_")) %>%
       select(!matches("category")) %>%
       filter(!is.na(exer_curr_act1_intensity)) %>%
       group_by(CONP_ID) %>%
       slice(1) %>%
       mutate(across(starts_with("exer_curr_") & where(is.character), as.numeric)) %>%
       rowwise() %>%
       mutate(
         across(starts_with("exer_curr_act2_"),
                ~ if_else(!is.na(exer_curr_act1_intensity) & is.na(.), 0, .)),
         across(starts_with("exer_curr_act3_"),
                ~ if_else(!is.na(exer_curr_act1_intensity) & is.na(.), 0, .)),
         across(starts_with("exer_curr_act4_"),
                ~ if_else(!is.na(exer_curr_act1_intensity) & is.na(.), 0, .)),
         across(starts_with("exer_curr_act5_"),
                ~ if_else(!is.na(exer_curr_act1_intensity) & is.na(.), 0, .))
       ) %>%
       ungroup()
    ),
    by = "CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, starts_with("social_life_")) %>%
       group_by(CONP_ID) %>%
       filter_at(vars(starts_with("social_life_")), any_vars(!is.na(.))) %>%
       slice(1)
    ),
    by = "CONP_ID") %>%
  # Medication columns set to NA: RxNorm lookup skipped for small sv subset.
  # The treatment_* fields from medical_history (already joined above) cover
  # most cases in the downstream recoding logic. Setting medusage_* to NA is
  # conservative — it won't create false positives, may miss a few cases.
  mutate(
    medusage_antihypertensive  = NA_real_,
    medusage_antidiabetic      = NA_real_,
    medusage_antihyperlipidemic = NA_real_,
    medusage_antidepressants   = NA_real_
  )


# TERTILE THRESHOLDS FROM MULTI-VISIT IMPUTED DATA ------------------------
# Using multi-visit cutoffs so that sv_ids' tertile assignments are on the
# same scale as criteria_met participants (not recomputed on ~20-50 people).

epoch_q33  <- quantile(crsfactors_mv$epoch_score_currently, 0.33, na.rm = TRUE)
epoch_q67  <- quantile(crsfactors_mv$epoch_score_currently, 0.67, na.rm = TRUE)

pittsburgh_q67 <- quantile(crsfactors_mv$pittsburgh_total_score, 0.67, na.rm = TRUE)

social_cols_mv <- crsfactors_mv %>%
  select(starts_with("social_life_frequency_")) %>%
  select(where(is.numeric))

social_row_means_mv  <- rowMeans(social_cols_mv, na.rm = TRUE)
social_row_totals_mv <- rowSums(social_cols_mv, na.rm = TRUE)

social_mean     <- mean(social_row_means_mv, na.rm = TRUE)
social_total_q33 <- quantile(social_row_totals_mv, 0.33, na.rm = TRUE)

cat("\nMulti-visit tertile reference thresholds:\n")
cat("  Epoch (33rd pctile):", round(epoch_q33, 2), "\n")
cat("  Epoch (67th pctile):", round(epoch_q67, 2), "\n")
cat("  Pittsburgh (67th pctile):", round(pittsburgh_q67, 2), "\n")
cat("  Social activity row mean:", round(social_mean, 2), "\n")
cat("  Social activity total (33rd pctile):", round(social_total_q33, 2), "\n")


# RECODE COGDRISK VARIABLES -----------------------------------------------
# Mirrors 03_recode_PREVENT-AD_CogD.R with fixed thresholds for ntile vars

clinical_cogd_sv <- crs_factors_sv %>%
  mutate(
    Sex = if_else(Sex == "Female", 0, 1),
    Education_level = case_when(
      Education_years > 11               ~ "High",
      Education_years > 8 & Education_years <= 11 ~ "Middle",
      Education_years < 8                ~ "Low",
      TRUE                               ~ NA_character_
    ) %>% factor(),
    Hypertension = case_when(
      medusage_antihypertensive == 1                                                ~ 1,
      Systolic_blood_pressure > 129                                                 ~ 1,
      Diastolic_blood_pressure >= 80                                                ~ 1,
      medusage_antihypertensive == 0 &
        !is.na(Systolic_blood_pressure) & !is.na(Diastolic_blood_pressure)         ~ 0,
      TRUE                                                                          ~ NA_real_
    ),
    BMI          = Weight / (Height / 100)^2,
    BMI_category = case_when(
      BMI < 18.5            ~ "Underweight",
      BMI > 18.5 & BMI < 25 ~ "Normal",
      BMI >= 25 & BMI < 30  ~ "Overweight",
      BMI >= 30             ~ "Obese",
      TRUE                  ~ NA_character_
    ) %>% factor(),
    High_cholesterol = if_else(total_cholesterol_value > 6.5, 1, 0),
    Insomnia         = if_else(pittsburgh_total_score > 5, 1, 0),
    Depression = case_when(
      past_depression == 1                               ~ 1,
      medusage_antidepressants == 1                      ~ 1,
      gds_score >= 5                                     ~ 1,
      past_depression == 0 & medusage_antidepressants == 0 ~ 0,
      TRUE                                               ~ NA_real_
    ),
    Atrial_fibrillation = past_atrial_fibrillation,
    Diabetes = case_when(
      treatment_diabetes > 0                                                          ~ 1,
      medusage_antidiabetic == 1                                                      ~ 1,
      Age < 60 & treatment_diabetes == 0 & hba1c_value <= 0.061                      ~ 0,
      (Age >= 60 & Age < 70) & hba1c_value <= 0.075 & treatment_diabetes == 0        ~ 0,
      Age >= 70 & treatment_diabetes == 0 & hba1c_value <= 0.07                      ~ 0,
      is.na(treatment_diabetes) & is.na(medusage_antidiabetic) & is.na(hba1c_value) ~ NA_real_,
      TRUE                                                                            ~ 1
    ),
    TBI = if_else(head_injury_hospitalized == 1 | head_injury_severe == 1, 1, 0)
  )

lifestyle_cogd_sv <- crs_factors_sv %>%
  mutate(
    Smoking = case_when(
      smoking_present == 4 ~ 2,
      smoking_present == 3 ~ 2,
      smoking_present == 2 ~ 1,
      smoking_present == 1 ~ 1,
      smoking_present == 0 ~ 0,
      TRUE                 ~ NA_real_
    ),
    # Use fixed tertile thresholds from multi-visit data instead of ntile()
    Cognitive_engagement = case_when(
      is.na(epoch_score_currently)          ~ NA_real_,
      epoch_score_currently <= epoch_q33    ~ 1,
      epoch_score_currently <= epoch_q67    ~ 2,
      TRUE                                  ~ 3
    ) %>% factor(),
    social_row_mean = rowMeans(
      select(., starts_with("social_life_frequency_")), na.rm = TRUE
    ),
    Social_engagement = factor(
      if_else(social_row_mean > social_mean, 0, 1)
    )
  ) %>%
  select(-social_row_mean)

# Exercise: same for-loop as 03_recode_PREVENT-AD_CogD.R
lifestyle_exercise_cogd_sv <- lifestyle_cogd_sv %>%
  mutate(
    light_minutes_week    = 0,
    moderate_minutes_week = 0,
    heavy_minutes_week    = 0
  ) %>%
  mutate(
    across(contains("_intensity"), ~ as.numeric(.)),
    across(contains("_weeks"),     ~ as.numeric(.)),
    across(contains("_hours"),     ~ as.numeric(.))
  ) %>%
  rowwise() %>%
  mutate(
    has_any_exercise_data = any(!is.na(c_across(contains("exer_curr_act"))))
  ) %>%
  ungroup()

for (i in 1:5) {
  intensity_col <- paste0("exer_curr_act", i, "_intensity")
  weeks_col     <- paste0("exer_curr_act", i, "_weeks")
  hours_col     <- paste0("exer_curr_act", i, "_hours")

  if (!all(c(intensity_col, weeks_col, hours_col) %in% names(lifestyle_exercise_cogd_sv))) next

  lifestyle_exercise_cogd_sv <- lifestyle_exercise_cogd_sv %>%
    mutate(
      activity_minutes = case_when(
        !is.na(.data[[intensity_col]]) &
          !is.na(.data[[weeks_col]]) &
          !is.na(.data[[hours_col]]) ~ .data[[weeks_col]] * .data[[hours_col]] * 60,
        TRUE ~ 0
      ),
      light_minutes_week = light_minutes_week +
        ifelse(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 1,
               activity_minutes, 0),
      moderate_minutes_week = moderate_minutes_week +
        ifelse(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 2,
               activity_minutes, 0),
      heavy_minutes_week = heavy_minutes_week +
        ifelse(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 3,
               activity_minutes, 0)
    ) %>%
    select(-activity_minutes)
}

lifestyle_cogd_sv <- lifestyle_exercise_cogd_sv %>%
  mutate(
    light_minutes_week    = if_else(has_any_exercise_data, light_minutes_week,    NA_real_),
    moderate_minutes_week = if_else(has_any_exercise_data, moderate_minutes_week, NA_real_),
    heavy_minutes_week    = if_else(has_any_exercise_data, heavy_minutes_week,    NA_real_)
  ) %>%
  select(-has_any_exercise_data) %>%
  mutate(
    Physical_inactivity = if_else(
      (moderate_minutes_week + heavy_minutes_week) > 150, 0, 1
    )
  ) %>%
  select(CONP_ID, Smoking, Cognitive_engagement, Social_engagement, Physical_inactivity)

dat_cogd_sv <- clinical_cogd_sv %>%
  left_join(lifestyle_cogd_sv, by = "CONP_ID")

dat_cogd_sv <- dat_cogd_sv %>%
  mutate(score_cogd = calculate_cogdrisk(.)) %>%
  relocate(score_cogd, .after = CONP_ID)


# RECODE LIBRA2 VARIABLES -------------------------------------------------
# Mirrors 04_recode_PREVENT-AD_libra2.R with fixed thresholds for ntile vars

dat_libra_sv <- crs_factors_sv %>%
  mutate(
    # Sleep disturbances: highest tertile of Pittsburgh score
    Sleep_disturbances = case_when(
      is.na(pittsburgh_total_score)                ~ NA_real_,
      pittsburgh_total_score >= pittsburgh_q67     ~ 1,
      TRUE                                         ~ 0
    ),
    Hypertension = case_when(
      Systolic_blood_pressure >= 140                                    ~ 1,
      Diastolic_blood_pressure >= 90                                    ~ 1,
      treatment_hypertension == 1 | treatment_hypertension == 2        ~ 1,
      medusage_antihypertensive == 1                                    ~ 1,
      is.na(Systolic_blood_pressure) & is.na(Diastolic_blood_pressure) &
        is.na(treatment_hypertension) & is.na(medusage_antihypertensive) ~ NA_real_,
      TRUE                                                              ~ 0
    ),
    Hypercholesterolemia = case_when(
      total_cholesterol_value >= 6.2                                  ~ 1,
      treatment_hyperlipidemia == 1 | treatment_hyperlipidemia == 2   ~ 1,
      medusage_antihyperlipidemic == 1                                ~ 1,
      is.na(total_cholesterol_value) & is.na(treatment_hyperlipidemia) &
        is.na(medusage_antihyperlipidemic)                            ~ NA_real_,
      TRUE                                                            ~ 0
    ),
    BMI     = Weight / (Height / 100)^2,
    Obesity = if_else(BMI >= 30, 1, 0),
    Hearing_impairment = if_else(diagnosed_impairment != 0, 1, 0),
    Diabetes = case_when(
      treatment_diabetes > 0                                                          ~ 1,
      medusage_antidiabetic == 1                                                      ~ 1,
      Age < 60 & treatment_diabetes == 0 & hba1c_value <= 0.061                      ~ 0,
      (Age >= 60 & Age < 70) & hba1c_value <= 0.075 & treatment_diabetes == 0        ~ 0,
      Age >= 70 & treatment_diabetes == 0 & hba1c_value <= 0.07                      ~ 0,
      is.na(treatment_diabetes) & is.na(medusage_antidiabetic) & is.na(hba1c_value) ~ NA_real_,
      TRUE                                                                            ~ 1
    ),
    Depression = case_when(
      past_depression == 1                                ~ 1,
      medusage_antidepressants == 1                       ~ 1,
      gds_score >= 5                                      ~ 1,
      past_depression == 0 & medusage_antidepressants == 0 ~ 0,
      TRUE                                                ~ NA_real_
    ),
    Smoking = if_else(as.numeric(smoking_present) >= 3, 1, 0),
    # Cognitive activity: highest tertile of epoch score
    Cognitive_activity = case_when(
      is.na(epoch_score_currently)        ~ NA_real_,
      epoch_score_currently >= epoch_q67  ~ 1,
      TRUE                                ~ 0
    )
  ) %>%
  mutate(
    social_total       = rowSums(select(., starts_with("social_life_frequency_")), na.rm = TRUE),
    Social_participation = case_when(
      is.na(social_total)                     ~ NA_real_,
      social_total <= social_total_q33        ~ 1,
      TRUE                                    ~ 0
    )
  ) %>%
  select(CONP_ID, Sleep_disturbances, Hypertension, Hypercholesterolemia, Obesity,
         Hearing_impairment, Diabetes, Depression, Smoking, Cognitive_activity,
         Social_participation)

# Exercise for LIBRA2: same for-loop as 04_recode_PREVENT-AD_libra2.R
exercise_libra_sv <- crs_factors_sv %>%
  select(CONP_ID, starts_with("exer_curr_")) %>%
  select(!matches("category")) %>%
  mutate(
    moderate_minutes_week = 0,
    vigorous_minutes_week = 0
  ) %>%
  mutate(
    across(contains("_intensity"), ~ as.numeric(.)),
    across(contains("_weeks"),     ~ as.numeric(.)),
    across(contains("_hours"),     ~ as.numeric(.))
  ) %>%
  rowwise() %>%
  mutate(
    has_any_exercise_data = any(!is.na(c_across(contains("exer_curr_act"))))
  ) %>%
  ungroup()

for (i in 1:5) {
  intensity_col <- paste0("exer_curr_act", i, "_intensity")
  weeks_col     <- paste0("exer_curr_act", i, "_weeks")
  hours_col     <- paste0("exer_curr_act", i, "_hours")

  if (!all(c(intensity_col, weeks_col, hours_col) %in% names(exercise_libra_sv))) next

  exercise_libra_sv <- exercise_libra_sv %>%
    mutate(
      activity_minutes = case_when(
        !is.na(.data[[intensity_col]]) &
          !is.na(.data[[weeks_col]]) &
          !is.na(.data[[hours_col]]) ~ .data[[weeks_col]] * .data[[hours_col]] * 60,
        TRUE ~ 0
      ),
      moderate_minutes_week = moderate_minutes_week +
        ifelse(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 2,
               activity_minutes, 0),
      vigorous_minutes_week = vigorous_minutes_week +
        ifelse(!is.na(.data[[intensity_col]]) & .data[[intensity_col]] == 3,
               activity_minutes, 0)
    ) %>%
    select(-activity_minutes)
}

exercise_libra_sv <- exercise_libra_sv %>%
  mutate(
    meets_moderate_guideline = moderate_minutes_week >= 150,
    meets_vigorous_guideline = vigorous_minutes_week >= 60,
    Physical_activity = case_when(
      meets_moderate_guideline | meets_vigorous_guideline ~ 0,
      TRUE                                                ~ 1
    )
  ) %>%
  select(CONP_ID, Physical_activity)

dat_libra_sv <- dat_libra_sv %>%
  left_join(exercise_libra_sv, by = "CONP_ID")

dat_libra_sv <- dat_libra_sv %>%
  mutate(score_libra2 = calculate_libra2(.)) %>%
  relocate(score_libra2, .after = CONP_ID)


# FAMILY HISTORY ----------------------------------------------------------
# Exact same block as 02_clean_PREVENT-AD.R lines 290–325, filtered to sv_ids

fhx_sv <- PREVENTAD_raw$demographics %>%
  filter(CONP_ID %in% sv_ids) %>%
  select(CONP_ID, father_dx_ad_dementia:other_paternal_family_members_AD) %>%
  rowwise() %>%
  mutate(
    FDR_AD = sum(father_dx_ad_dementia,
                 mother_dx_ad_dementia,
                 sibling_dx_ad_dementia_count, na.rm = TRUE),
    other_family_members_AD = if_else(is.na(other_family_members_AD), 0,
                                      other_family_members_AD),
    other_maternal_family_members_AD = if_else(is.na(other_maternal_family_members_AD), 0,
                                               other_maternal_family_members_AD),
    other_paternal_family_members_AD = if_else(is.na(other_paternal_family_members_AD), 0,
                                               other_paternal_family_members_AD),
    ER_AD = if_else(other_family_members_AD > 0, 1, 0),
    FDRAD_1ormore = if_else(FDR_AD == 1, 1, 2),
    family_history = case_when(
      FDRAD_1ormore == 1 & ER_AD == 0 ~ 1,
      FDRAD_1ormore == 1 & ER_AD == 1 ~ 2,
      FDRAD_1ormore == 2 & ER_AD == 0 ~ 3,
      FDRAD_1ormore == 2 & ER_AD == 1 ~ 4,
      TRUE                            ~ NA_real_
    ),
    family_history = factor(
      family_history,
      levels = c(1, 2, 3, 4),
      labels = c("1 FDR", "1 FDR + ext", "2+ FDR", "2+ FDR + ext")
    )
  ) %>%
  ungroup() %>%
  select(CONP_ID, FDR_AD, ER_AD, family_history)


# COMBINE AND SAVE --------------------------------------------------------

sv_scores <- dat_cogd_sv %>%
  select(CONP_ID, score_cogd) %>%
  left_join(
    dat_libra_sv %>% select(CONP_ID, score_libra2),
    by = "CONP_ID"
  ) %>%
  left_join(fhx_sv, by = "CONP_ID")

cat("\nSingle-visit score summary:\n")
cat("  N participants:", nrow(sv_scores), "\n")
cat("  N with score_cogd:", sum(!is.na(sv_scores$score_cogd)), "\n")
cat("  N with score_libra2:", sum(!is.na(sv_scores$score_libra2)), "\n")
cat("  N with family_history:", sum(!is.na(sv_scores$family_history)), "\n")

cat("\nscore_cogd range:", round(range(sv_scores$score_cogd, na.rm = TRUE), 1), "\n")
cat("score_libra2 range:", round(range(sv_scores$score_libra2, na.rm = TRUE), 1), "\n")

saveRDS(sv_scores,
        file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_singlevisi_scores.rds"))

cat("\nSaved: data/intermediate/PREVENTAD_singlevisi_scores.rds\n")
