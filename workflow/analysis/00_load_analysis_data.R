
# SETUP -------------------------------------------------------------------

# Load files
source('workflow/scripts/config.R')

apoe_dat <- readRDS(file.path(DATA_CLEANED_PATH$filtered, "PREVENTAD_apoe_filtered.rds"))
biomarkers <- readRDS(file.path(DATA_CLEANED_PATH$filtered, 
                                "PREVENTAD_biomarkers_filtered.rds"))
clinical <- readRDS(file.path(DATA_CLEANED_PATH$filtered, "PREVENTAD_clinical_imp_filtered_raw.rds"))
cogd_dat <- readRDS(file.path(DATA_CLEANED_PATH$cleaned, "PREVENTAD_cogd_scored_dat.rds"))
diagnosis <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_diagnosis.rds"))
fhx <- readRDS(file.path(DATA_CLEANED_PATH$filtered, 
                         "PREVENTAD_fhx_filtered.rds"))
grs <- readRDS(file.path(DATA_CLEANED_PATH$cleaned, "PREVENTAD_GRS.rds"))
libra_dat <- readRDS(file.path(DATA_CLEANED_PATH$cleaned, "PREVENTAD_libra_scored_dat.rds"))
PREVENTAD_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_dat.rds"))


# PREDICTORS DATAFRAME ------------------------------------------------------
# Create predictors dataframe
predictors_dat <- grs %>%
  left_join((clinical %>%
               select(CONP_ID, Candidate_Age, Sex, Weight, Height) %>%
               mutate(Candidate_Age = Candidate_Age/12,
                      BMI = Weight/(Height/100)^2) %>%
               rename(age = Candidate_Age,
                      sex = Sex) %>%
               select(-Weight, -Height)),
            by="CONP_ID") %>%
  left_join(apoe_dat, by="CONP_ID") %>%
  left_join((cogd_dat %>% 
               select(CONP_ID, score_cogd)),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$demographics %>% 
               select(CONP_ID, educ = Education_years)),
            by="CONP_ID")
  left_join(fhx, by="CONP_ID") %>%
  left_join((libra_dat %>% 
               select(CONP_ID, score_libra2)),
            by="CONP_ID") %>%
  select(CONP_ID, age, sex, BMI, apoe, 
         score_cogd, educ, GRS, score_libra2, family_history)

# Save processed dataframe
saveRDS(predictors_dat, file.path(DATA_CLEANED_PATH$filtered, "PREVENTAD_predictors.rds"))

# BIOMARKER DATAFRAME -----------------------------------------------------
# Baseline biomarker dataframe
baseline_dat <- biomarkers %>%
  group_by(CONP_ID) %>%
  select(CONP_ID, all_of(BIOMARKERS)) %>%
  slice_head() %>%
  left_join(predictors_dat, by = "CONP_ID")

# Biomarker dataframe (last visit) 
ptau217_last <- biomarkers %>%
  group_by(CONP_ID) %>%
  filter(!is.na(ptau217)) %>%
  select(CONP_ID, Date_taken, age, ptau217) %>%
  slice_tail() 

other_biomarkers_last <- biomarkers %>%
  group_by(CONP_ID) %>%
  select(CONP_ID, all_of(BIOMARKERS)) %>%
  select(-ptau217) %>%
  filter_at(vars(-CONP_ID), all_vars(!is.na(.))) %>%
  slice_tail() 

biomarker_last_dat <- predictors_dat %>% 
  right_join((ptau217_last %>%
               select(CONP_ID, ptau217)), by="CONP_ID") %>%
  left_join(other_biomarkers_last, by="CONP_ID") %>%
  mutate(gfap_log = log(gfap),
         nfl_log = log(nfl),
         ptau181_log = log(ptau181),
         ptau217_log = log(ptau217),
         ptau231_sqrt = sqrt(ptau231),
         ptau217_ab42_ratio_log = log(ptau217_ab42_ratio)) 

# Longitudinal biomarker dataframe
all_dat <- biomarkers %>%
  left_join(predictors_dat, by = "CONP_ID") %>%
  mutate(gfap_log = log(gfap),
         nfl_log = log(nfl),
         ptau181_log = log(ptau181),
         ptau217_log = log(ptau217),
         ptau231_sqrt = sqrt(ptau231),
         ptau217_ab42_ratio_log = log(ptau217_ab42_ratio)) 

# Save processed dataframes
saveRDS(baseline_dat, file.path(DATA_CLEANED_PATH$filtered, "PREVENTAD_baseline.rds"))
saveRDS(biomarker_last_dat, file.path(DATA_CLEANED_PATH$filtered, "PREVENTAD_lastvisit.rds"))
saveRDS(all_dat, file.path(DATA_CLEANED_PATH$filtered, "PREVENTAD_longitudinal.rds"))
