
# SETUP -------------------------------------------------------------------

# Load files
source('workflow/scripts/config.R')
source(CRS_FN)

BRANCH_cleaned_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH$BRANCH, "BRANCH_cleaned.rds"))


# RECODING ----------------------------------------------------------------

# Recode clinical data
clinical_cogd <- BRANCH_cleaned_dat$demographics %>%
  left_join(
    BRANCH_cleaned_dat$physical %>% select(PIDN, height, weight, bpsys, bpdias),
    by = "PIDN"
  ) %>%
  left_join(
    BRANCH_cleaned_dat$health_history %>% select(PIDN, diabetes, dep2yrs, quitsmok, smokyrs,
                                                 tbi, cvafib, insomn),
    by = "PIDN"
  ) %>%
  left_join(
    BRANCH_cleaned_dat$clinical_labs %>% select(PIDN, total_cholesterol, hba1c),
    by = "PIDN"
  ) %>%
  left_join(
    BRANCH_cleaned_dat$APOE %>% select(PIDN, apoe),
    by = "PIDN"
  ) %>%
  mutate(
    Sex = if_else(gender == "Female", 0, 1) %>% factor(), # female == 0, male == 1
    Age = age_at_DCDate,
    Education_level = case_when(
      educ > 11 ~ "High",
      educ > 8 & educ <= 11 ~ "Middle",
      educ <= 8 ~ "Low",
      TRUE ~ NA) %>% factor(),
    
    # Hypertension risk factor present == 1
    Hypertension = case_when(
      bpsys > 129 ~ 1,
      bpdias >= 80 ~ 1,
      !is.na(bpsys) & !is.na(bpdias) ~ 0,
      TRUE ~ NA) %>% factor(),
    
    # Calculate BMI (height in inches, weight in pounds)
    BMI = (weight / height^2) * 703,
    BMI_category = case_when(
      BMI < 18.5 ~ "Underweight",
      BMI >= 18.5 & BMI < 25 ~ "Normal",
      BMI >= 25 & BMI < 30 ~ "Overweight",
      BMI >= 30 ~ "Obese",
      TRUE ~ NA) %>% factor(),
    
    # High cholesterol risk factor present == 1
    High_cholesterol = if_else(total_cholesterol > 250, 1, 0) %>% factor(),
    
    # Insomnia risk factor present == 1
    Insomnia = if_else(insomn == 1, 1, 0) %>% factor(),
    
    # Depression risk factor present == 1
    Depression = if_else(dep2yrs == 1, 1, 0) %>% factor(),
    
    # Atrial fibrillation
    Atrial_fibrillation = if_else(cvafib == 1, 1, 0) %>% factor(),
    
    # Diabetes risk factor present == 1
    Diabetes = case_when(
      diabetes == 1 ~ 1,
      Age < 60 & diabetes == 0 & hba1c <= 6.1 ~ 0,
      (Age >= 60 & Age < 70) & hba1c <= 7.5 & diabetes == 0 ~ 0,
      Age >= 70 & diabetes == 0 & hba1c <= 7.0 ~ 0,
      is.na(diabetes) & is.na(hba1c) ~ NA,
      TRUE ~ 1) %>% factor(),
    
    # TBI
    TBI = case_when(
      tbi == 2 ~ 1,  # Yes (multiple)
      tbi == 1 ~ 1,  # Yes
      tbi == 0 ~ 0,  # No
      TRUE ~ NA) %>% factor()
  ) %>%
  select(PIDN, Sex, Age, Education_level, Hypertension, BMI, BMI_category,
         High_cholesterol, Insomnia, Depression, Atrial_fibrillation,
         Diabetes, TBI, apoe)

# Recode lifestyle data
lifestyle_cogd <- BRANCH_cleaned_dat$health_history %>%
  left_join(BRANCH_cleaned_dat$pase, by = "PIDN") %>%
  left_join(BRANCH_cleaned_dat$sni, by = "PIDN") %>%
  left_join(BRANCH_cleaned_dat$diet, by = "PIDN") %>%
  left_join(BRANCH_cleaned_dat$cas, by = "PIDN") %>%
  mutate(
    # Smoking: 2 == current, 1 == former, 0 == non smoker
    Smoking = case_when(
      quitsmok == 0 & smokyrs > 0 ~ 2,  # Never quit but has smoking years = current
      quitsmok == 1 ~ 1,  # Former smoker
      is.na(smokyrs) | smokyrs == 0 ~ 0,  # Non-smoker
      TRUE ~ NA),
    
    # Diet: MIND score tertiles (higher tertile = better diet = lower risk)
    # For CogDRisk, higher tertile should be protective
    Diet = ntile(mind_score, 3) %>% factor(),
    
    # Cognitive engagement: CAS total points tertiles (higher tertile = more engagement)
    Cognitive_engagement = ntile(cas_tot_pts, 3) %>% factor(),
    
    # Social engagement based on SNI
    # 0 == not lonely (high social contact), 1 == lonely (low social contact)
    Social_engagement = case_when(
      SNI_HighContact >= 4 ~ 0,  # High contact = not lonely
      SNI_HighContact < 4 ~ 1,   # Low contact = lonely
      TRUE ~ NA) %>% factor(),
    
    # Physical inactivity based on PASE
    # Using median split: below median = inactive (1), above median = active (0)
    Physical_inactivity = if_else(
      pase_total < median(pase_total, na.rm = TRUE), 1, 0
    ) %>% factor()
  ) %>%
  select(PIDN, Smoking, Diet, Cognitive_engagement, Social_engagement, Physical_inactivity)

# Merge
dat_cogd <- clinical_cogd %>%
  left_join(lifestyle_cogd, by = "PIDN") %>%
  ungroup()

# Scoring
dat_cogd_scored <- dat_cogd %>%
  mutate(score_cogd = calculate_cogdrisk(.)) %>%
  relocate(score_cogd, .after = PIDN)

# Save
saveRDS(dat_cogd,
        file.path(DATA_INTERMEDIATE_PATH$BRANCH, "BRANCH_cogd_dat.rds"))
saveRDS(dat_cogd_scored,
        file.path(DATA_INTERMEDIATE_PATH$BRANCH, "BRANCH_cogd_scored_dat.rds"))


# IMPUTED DATA ------------------------------------------------------------

# Use imputed baseline data
baseline_imp <- BRANCH_cleaned_dat$baseline_imputed

# Recode clinical data with imputed values
clinical_cogd_imp <- baseline_imp %>%
  select(PIDN, age_at_DCDate, gender, educ, height, weight, bpsys, bpdias,
         diabetes, dep2yrs, quitsmok, smokyrs, tbi, cvafib, insomn,
         total_cholesterol, hba1c) %>%
  left_join(BRANCH_cleaned_dat$APOE %>% select(PIDN, apoe), by = "PIDN") %>%
  mutate(
    Sex = if_else(gender == "Female", 0, 1) %>% factor(),
    Age = age_at_DCDate,
    Education_level = case_when(
      educ > 11 ~ "High",
      educ > 8 & educ <= 11 ~ "Middle",
      educ <= 8 ~ "Low",
      TRUE ~ NA) %>% factor(),
    Hypertension = case_when(
      bpsys > 129 ~ 1,
      bpdias >= 80 ~ 1,
      !is.na(bpsys) & !is.na(bpdias) ~ 0,
      TRUE ~ NA) %>% factor(),
    BMI = (weight / height^2) * 703,
    BMI_category = case_when(
      BMI < 18.5 ~ "Underweight",
      BMI >= 18.5 & BMI < 25 ~ "Normal",
      BMI >= 25 & BMI < 30 ~ "Overweight",
      BMI >= 30 ~ "Obese",
      TRUE ~ NA) %>% factor(),
    High_cholesterol = if_else(total_cholesterol > 250, 1, 0) %>% factor(),
    Insomnia = if_else(insomn == 1, 1, 0) %>% factor(),
    Depression = if_else(dep2yrs == 1, 1, 0) %>% factor(),
    Atrial_fibrillation = if_else(cvafib == 1, 1, 0) %>% factor(),
    Diabetes = case_when(
      diabetes == 1 ~ 1,
      Age < 60 & diabetes == 0 & hba1c <= 6.1 ~ 0,
      (Age >= 60 & Age < 70) & hba1c <= 7.5 & diabetes == 0 ~ 0,
      Age >= 70 & diabetes == 0 & hba1c <= 7.0 ~ 0,
      is.na(diabetes) & is.na(hba1c) ~ NA,
      TRUE ~ 1) %>% factor(),
    TBI = case_when(
      tbi == 2 ~ 1,
      tbi == 1 ~ 1,
      tbi == 0 ~ 0,
      TRUE ~ NA) %>% factor()
  ) %>%
  select(PIDN, Sex, Age, Education_level, Hypertension, BMI, BMI_category,
         High_cholesterol, Insomnia, Depression, Atrial_fibrillation,
         Diabetes, TBI, apoe)

# Recode lifestyle data with imputed values
lifestyle_cogd_imp <- baseline_imp %>%
  select(PIDN, quitsmok, smokyrs, pase_total, SNI_HighContact, mind_score, cas_tot_pts) %>%
  mutate(
    Smoking = case_when(
      quitsmok == 0 & smokyrs > 0 ~ 2,
      quitsmok == 1 ~ 1,
      is.na(smokyrs) | smokyrs == 0 ~ 0,
      TRUE ~ NA),
    Diet = ntile(mind_score, 3) %>% factor(),
    Cognitive_engagement = ntile(cas_tot_pts, 3) %>% factor(),
    Social_engagement = case_when(
      SNI_HighContact >= 4 ~ 0,
      SNI_HighContact < 4 ~ 1,
      TRUE ~ NA) %>% factor(),
    Physical_inactivity = if_else(
      pase_total < median(pase_total, na.rm = TRUE), 1, 0
    ) %>% factor()
  ) %>%
  select(PIDN, Smoking, Diet, Cognitive_engagement, Social_engagement, Physical_inactivity)

# Merge imputed data
dat_cogd_imp <- clinical_cogd_imp %>%
  left_join(lifestyle_cogd_imp, by = "PIDN") %>%
  ungroup()

# Score imputed data
dat_cogd_imp_scored <- dat_cogd_imp %>%
  mutate(score_cogd = calculate_cogdrisk(.)) %>%
  relocate(score_cogd, .after = PIDN)

# Save imputed datasets
saveRDS(dat_cogd_imp,
        file.path(DATA_INTERMEDIATE_PATH$BRANCH, "BRANCH_cogd_imp_dat.rds"))
saveRDS(dat_cogd_imp_scored,
        file.path(DATA_INTERMEDIATE_PATH$BRANCH, "BRANCH_cogd_imp_scored_dat.rds"))
