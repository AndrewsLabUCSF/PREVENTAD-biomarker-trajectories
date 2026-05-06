
# SETUP -------------------------------------------------------------------

# Load files
source('workflow/scripts/config.R')

# Load intermediate datasets from data wrangling pipeline
biomarkers <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_biomarkers.rds"))
apoe_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_APOE.rds"))
fhx_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_fhx.rds"))
grs_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_GRS.rds"))
mci_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_MCI.rds"))

# Load CRS factors and scored risk scores
crs_factors <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_crsfactors_imputed.rds"))
cogd_scored <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_cogd_scored.rds"))
libra_scored <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_libra2_scored.rds"))

# IDs of participants with >25% missing CRS variables
missing_ids <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_missingcrs.rds"))


# PREDICTORS DATAFRAME ------------------------------------------------------
# Create predictors dataframe by combining all risk factor sources and filter out IDs
predictors_dat <- crs_factors %>%
  filter(!CONP_ID %in% missing_ids$CONP_ID) %>%
  select(CONP_ID, age=Age, sex=Sex, educ=Education_years) %>%
  left_join((grs_dat %>% select(CONP_ID, GRS)), by = "CONP_ID") %>%
  left_join((apoe_dat %>% select(CONP_ID, apoe)), by = "CONP_ID") %>%
  left_join((fhx_dat %>% select(CONP_ID, FDR_AD)), by = "CONP_ID") %>%
  left_join((cogd_scored %>% select(CONP_ID, score_cogd)), by = "CONP_ID") %>%
  left_join((libra_scored %>% select(CONP_ID, score_libra2)), by = "CONP_ID") %>%
  select(CONP_ID, age, sex, educ, apoe, score_cogd, GRS, score_libra2, FDR_AD) %>%
  # Risk burden components: 
  mutate(
    grs_quartile      = ntile(GRS, 4),
    grs_burden        = if_else(grs_quartile == 4, 1, 0),
    fhx_burden        = if_else(FDR_AD > 1, 1, 0),
    # fhx_burden        = if_else(family_history == "1 FDR", 0, 1),
    apoe_burden       = if_else(apoe == "e4+", 1, 0),
    
    # CRS cutoff strategy: 1 SD above mean (selected via aim1_crs_cutoffs analysis)
    cogdrisk_burden   = if_else(score_cogd > mean(score_cogd, na.rm = TRUE) +
                                  sd(score_cogd, na.rm = TRUE), 1, 0),
    libra2_burden     = if_else(score_libra2 > mean(score_libra2, na.rm = TRUE) +
                                  sd(score_libra2, na.rm = TRUE), 1, 0),
    risk_score_cogd   = rowSums(across(c(fhx_burden, apoe_burden,
                                         cogdrisk_burden, grs_burden)), na.rm = TRUE),
    risk_score_libra2 = rowSums(across(c(fhx_burden, apoe_burden,
                                         libra2_burden, grs_burden)), na.rm = TRUE)
  ) 

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

baseline_dat <- predictors_dat %>%
  left_join((ptau217_first %>% select(CONP_ID, ptau217)), by = "CONP_ID") %>%
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

lastvisit_dat <- predictors_dat %>%
  left_join((ptau217_last %>% select(CONP_ID, ptau217)), by = "CONP_ID") %>%
  left_join(other_biomarkers_last, by = "CONP_ID") %>%
  mutate(
    # Add log transformations
    gfap_log = log(gfap),
    nfl_log = log(nfl),
    ptau181_log = log(ptau181),
    ptau217_log = log(ptau217),
    ptau217_ab42_ratio_log = log(ptau217_ab42_ratio)
  )

# Filter out participants that progressed to MCI
lastvisit_nomci_dat <- lastvisit_dat %>% filter(!CONP_ID %in% mci_dat$CONP_ID)


## Longitudinal ----
mean_years <- mean(biomarkers$years, na.rm = TRUE)

all_dat <- predictors_dat %>%
  left_join(biomarkers, by = "CONP_ID") %>%
  mutate(
    # Add log transformations
    gfap_log = log(gfap),
    nfl_log = log(nfl),
    ptau181_log = log(ptau181),
    ptau217_log = log(ptau217),
    ptau217_ab42_ratio_log = log(ptau217_ab42_ratio),
    
    # Grand-mean centering time around the average follow-up time across all observations
    years_centered  = years - mean_years
  ) %>%
  select(-age.y) %>%
  rename(age = age.x)

# Filter out participants that progressed to MCI
all_nomci_dat <- all_dat %>% filter(!CONP_ID %in% mci_dat$CONP_ID)


# SAVE FINAL DATAFRAMES ---------------------------------------------------
saveRDS(baseline_dat, file.path(DATA_CLEANED_PATH, "PREVENTAD_baseline.rds"))
saveRDS(lastvisit_dat, file.path(DATA_CLEANED_PATH, "PREVENTAD_lastvisit_all.rds"))
saveRDS(lastvisit_nomci_dat, file.path(DATA_CLEANED_PATH, "PREVENTAD_lastvisit_nomci.rds"))
saveRDS(all_dat, file.path(DATA_CLEANED_PATH, "PREVENTAD_longitudinal.rds"))
saveRDS(all_nomci_dat, file.path(DATA_CLEANED_PATH, "PREVENTAD_longitudinal_nomci.rds"))

# Print summary
cat("\n=== Data Loading Complete ===\n")
cat("Baseline:", nrow(baseline_dat), "participants\n")
cat("Last visit:", nrow(lastvisit_dat), "participants\n")
cat("Longitudinal data:", nrow(all_dat), "observations from",
    length(unique(all_dat$CONP_ID)), "participants\n")

