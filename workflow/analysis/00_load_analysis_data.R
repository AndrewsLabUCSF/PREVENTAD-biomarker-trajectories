
# SETUP -------------------------------------------------------------------

# Load files
source('workflow/scripts/config.R')

# Load intermediate datasets from data wrangling pipeline
biomarkers <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_biomarkers.rds"))
apoe_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_APOE.rds"))
fhx_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_fhx.rds"))
grs_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_GRS.rds"))


# Load CRS factors and scored risk scores
crs_factors <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_crsfactors_imputed.rds"))
cogd_scored <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_cogd_scored.rds"))
libra_scored <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_libra2_scored.rds"))


# PREDICTORS DATAFRAME ------------------------------------------------------
# Create predictors dataframe by combining all risk factor sources
predictors_dat <- crs_factors %>%
  select(CONP_ID, age=Age, sex=Sex, educ=Education_years) %>%
  left_join((grs_dat %>% select(CONP_ID, GRS)), by = "CONP_ID") %>%
  left_join((apoe_dat %>% select(CONP_ID, apoe)), by = "CONP_ID") %>%
  left_join((fhx_dat %>% select(CONP_ID, family_history)), by = "CONP_ID") %>%
  left_join((cogd_scored %>% select(CONP_ID, score_cogd)), by = "CONP_ID") %>%
  left_join((libra_scored %>% select(CONP_ID, score_libra2)), by = "CONP_ID") %>%
  select(CONP_ID, age, sex, educ, apoe, score_cogd, GRS, score_libra2, family_history)

# Save predictors dataframe
saveRDS(predictors_dat, file.path(DATA_CLEANED_PATH, "PREVENTAD_predictors.rds"))


# BIOMARKER DATAFRAMES -----------------------------------------------------

## Baseline ----
# Handle ptau217 separately as it may have different availability
ptau217_first <- biomarkers %>%
  group_by(CONP_ID) %>%
  filter(!is.na(ptau217)) %>%
  select(CONP_ID, Date_taken, age, ptau217) %>%
  slice_head() %>%
  ungroup()

other_biomarkers_first <- biomarkers %>%
  group_by(CONP_ID) %>%
  select(CONP_ID, all_of(BIOMARKERS)) %>%
  select(-ptau217) %>%
  filter_at(vars(-CONP_ID), all_vars(!is.na(.))) %>%
  slice_head() %>%
  ungroup()

biomarker_baseline_dat <- predictors_dat %>%
  right_join((ptau217_first %>% select(CONP_ID, ptau217)), by = "CONP_ID") %>%
  left_join(other_biomarkers_first, by = "CONP_ID") %>%
  # Add log transformations
  mutate(
    gfap_log = log(gfap),
    nfl_log = log(nfl),
    ptau181_log = log(ptau181),
    ptau217_log = log(ptau217),
    ptau217_ab42_ratio_log = log(ptau217_ab42_ratio)
  )


## Last visit ----
# Handle ptau217 separately as it may have different availability
ptau217_last <- biomarkers %>%
  group_by(CONP_ID) %>%
  filter(!is.na(ptau217)) %>%
  select(CONP_ID, Date_taken, age, ptau217) %>%
  slice_tail() %>%
  ungroup()

other_biomarkers_last <- biomarkers %>%
  group_by(CONP_ID) %>%
  select(CONP_ID, all_of(BIOMARKERS)) %>%
  select(-ptau217) %>%
  filter_at(vars(-CONP_ID), all_vars(!is.na(.))) %>%
  slice_tail() %>%
  ungroup()

biomarker_last_dat <- predictors_dat %>%
  right_join((ptau217_last %>% select(CONP_ID, ptau217)), by = "CONP_ID") %>%
  left_join(other_biomarkers_last, by = "CONP_ID") %>%
  # Add log transformations
  mutate(
    gfap_log = log(gfap),
    nfl_log = log(nfl),
    ptau181_log = log(ptau181),
    ptau217_log = log(ptau217),
    ptau217_ab42_ratio_log = log(ptau217_ab42_ratio)
  )

## Longitudinal ----
all_dat <- biomarkers %>%
  left_join(predictors_dat, by = "CONP_ID") %>%
  # Add log transformations
  mutate(
    gfap_log = log(gfap),
    nfl_log = log(nfl),
    ptau181_log = log(ptau181),
    ptau217_log = log(ptau217),
    ptau217_ab42_ratio_log = log(ptau217_ab42_ratio)
  ) %>%
  select(-age.y) %>%
  rename(age = age.x)


# SAVE FINAL DATAFRAMES ---------------------------------------------------
saveRDS(biomarker_baseline_dat, file.path(DATA_CLEANED_PATH, "PREVENTAD_baseline.rds"))
saveRDS(biomarker_last_dat, file.path(DATA_CLEANED_PATH, "PREVENTAD_lastvisit.rds"))
saveRDS(all_dat, file.path(DATA_CLEANED_PATH, "PREVENTAD_longitudinal.rds"))

# Print summary
cat("\n=== Data Loading Complete ===\n")
cat("Baseline:", nrow(baseline_dat), "participants\n")
cat("Last visit:", nrow(biomarker_last_dat), "participants\n")
cat("Longitudinal data:", nrow(all_dat), "observations from",
    length(unique(all_dat$CONP_ID)), "participants\n")

