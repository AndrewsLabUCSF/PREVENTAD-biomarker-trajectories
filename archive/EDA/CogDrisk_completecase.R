
# SETUP -------------------------------------------------------------------
library(ggpubr)
library(missForest)

source('workflow/scripts/config.R')
source(CRS_FN)
source(COGD_RECODE_FN)

clinical_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_clinical_raw.rds"))
lifestyle_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_lifestyle_raw.rds"))
dat_cogd <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_cogd_dat.rds"))
dat_cogd_scored <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_cogd_scored.rds"))


# CALCULATING COGDRISK ----------------------------------------------------
data_cogd_scored <- data_cogd %>% 
  select(-c(Education_years, BMI, LDL_value)) %>%
  mutate(num_cogdrisk_vars = (rowSums(!is.na(.))) - 1) %>%
  mutate(complete_scores = complete.cases(.)) %>%
  mutate(score = calculate_cogdrisk(.)) %>%
  relocate(score, .after=CONP_ID) %>%
  relocate(complete_scores, .after=score) %>%
  relocate(num_cogdrisk_vars, .after=complete_scores) 

saveRDS(data_cogd_scored, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_cogd_scored.rds"))


# MISSING DATA IMPUTATION -------------------------------------------------

## Imputing all NA vars ---------------------------------------------------
clinical_vars_to_imp <- names(clinical_raw)[sapply(clinical_raw, anyNA)]
lifestyle_vars_to_imp <- names(lifestyle_raw)[sapply(lifestyle_raw, anyNA)]

clinical_imp <- clinical_raw %>%
  select(all_of(clinical_vars_to_imp)) %>%
  mutate_at(c("head_injury_hospitalized", "head_injury_severe"), as.factor)
lifestyle_imp <- lifestyle_raw %>%
  select(all_of(lifestyle_vars_to_imp)) %>%
  mutate(exer_curr_act2_days = as.numeric(exer_curr_act2_days)) %>%
  mutate(smoking_present = as.factor(smoking_present)) %>%
  mutate(across(ends_with("_intensity"), as.factor)) %>%
  mutate(across(starts_with("social_"), as.factor)) 

# Use missForest to impute NAs
clinical_imp_mf <- missForest(clinical_imp, verbose=TRUE, maxiter=10, ntree=100)
clinical_imp_mf_dat <- cbind(clinical_imp_mf$ximp, CONP_ID=clinical_raw$CONP_ID) %>%
  relocate(CONP_ID) 
lifestyle_imp_mf <- missForest(lifestyle_imp, verbose=TRUE, maxiter=10, ntree=100)
lifestyle_imp_mf_dat <- cbind(lifestyle_imp_mf$ximp, CONP_ID=lifestyle_raw$CONP_ID) %>%
  relocate(CONP_ID) %>%
  mutate(across(starts_with("social_"), as.numeric))

# Update dataframes with imputed data
clinical_imp_raw <- clinical_raw
clinical_imp_raw[names(clinical_imp_mf_dat)] <- clinical_imp_mf_dat
lifestyle_imp_raw <- lifestyle_raw
lifestyle_imp_raw[names(lifestyle_imp_mf_dat)] <- lifestyle_imp_mf_dat

# Recode imputed data to CogDrisk
cogd_dat_imp <- prepare_preventad_cogdrisk(clinical_imp_raw,
                                           lifestyle_imp_raw)

# Calculate CogD on recoded imputed data
cogd_dat_imp_scored <- cogd_dat_imp %>% 
  select(-c(Education_years, BMI, LDL_value)) %>%
  mutate(num_cogdrisk_vars = (rowSums(!is.na(.))) - 1) %>%
  mutate(complete_scores = complete.cases(.)) %>%
  mutate(score = calculate_cogdrisk(.)) %>%
  relocate(score, .after=CONP_ID) %>%
  relocate(complete_scores, .after=score) %>%
  relocate(num_cogdrisk_vars, .after=complete_scores) %>%
  select(-complete_scores, -num_cogdrisk_vars)

# Save intermediate data
saveRDS(cogd_dat_imp, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_cogd_imp_raw.rds"))
saveRDS(clinical_imp_raw, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_clinical_imp_raw.rds"))
saveRDS(lifestyle_imp_raw, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_lifestyle_imp_raw.rds"))
saveRDS(cogd_dat_imp_scored, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_cogd_imp_scored.rds"))


## Imputing only first exercise variables if missing ----------------------
# Selecting NA columns and first exercise activity vars
lifestyle_exer_imp <- lifestyle_raw %>%
  select(smoking_present, epoch_score_currently, starts_with("exer_curr_act1_"),
         starts_with("social_")) %>%
  mutate(smoking_present = as.factor(smoking_present)) %>%
  mutate(across(ends_with("_intensity"), as.factor)) %>%
  mutate(across(starts_with("social_"), as.factor)) 

# Use missForest to impute NAs
lifestyle_exer_imp_mf <- missForest(lifestyle_exer_imp, verbose=TRUE, maxiter=10, ntree=100)
lifestyle_exer_imp_mf_dat <- cbind(lifestyle_exer_imp_mf$ximp, CONP_ID=lifestyle_raw$CONP_ID) %>%
  relocate(CONP_ID) %>%
  mutate(across(starts_with("social_"), as.numeric))

# Update dataframe with imputed data
lifestyle_exer_imp_raw <- lifestyle_raw
lifestyle_exer_imp_raw[names(lifestyle_exer_imp_mf_dat)] <- lifestyle_exer_imp_mf_dat

# Recode imputed data to CogDrisk
cogd_exer_dat_imp <- prepare_preventad_cogdrisk(clinical_imp_raw,
                                                lifestyle_exer_imp_raw)

# Calculcate CogD on recoded imputed data
cogd_exer_dat_imp_scored <- cogd_exer_dat_imp %>% 
  select(-c(Education_years, BMI, LDL_value)) %>%
  mutate(num_cogdrisk_vars = (rowSums(!is.na(.))) - 1) %>%
  mutate(complete_scores = complete.cases(.)) %>%
  mutate(score = calculate_cogdrisk(.)) %>%
  relocate(score, .after=CONP_ID) %>%
  relocate(complete_scores, .after=score) %>%
  relocate(num_cogdrisk_vars, .after=complete_scores) %>%
  select(-complete_scores, -num_cogdrisk_vars)

# Save intermediate data
saveRDS(cogd_exer_dat_imp, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_cogd_exercise1_imp_raw.rds"))
saveRDS(cogd_exer_dat_imp_scored, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_cogd_exercise1_imp_scored.rds"))
