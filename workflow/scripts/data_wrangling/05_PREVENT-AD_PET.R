
# =============================================================================
# 05_PREVENT-AD_PET.R
# Process PREVENT-AD PET imaging data for amyloid positivity analysis
# =============================================================================

# SETUP -----------------------------------------------------------------------

source('workflow/scripts/config.R')

library(tidyverse)
library(lubridate)

# INCLUSION CRITERIA ----------------------------------------------------------
# Use pre-processed biomarker data (already filtered by inclusion criteria from 02_clean):
#   - CONP_ID with GWAS data
#   - CONP_ID with at least 2 biomarker visits

biomarkers <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_biomarkers.rds"))
criteria_met <- biomarkers %>% distinct(CONP_ID)

cat("=== INCLUSION CRITERIA (from pre-processed biomarker data) ===\n")
cat("Participants meeting criteria:", nrow(criteria_met), "\n\n")

# LOAD PET DATA ---------------------------------------------------------------

# NAV amyloid PET with whole cerebellum reference
pet_nav <- read_csv(
  file.path(PREVENTAD_DATA_PATH,
            "PET_NAV_SUVR_ref-wholeCerebellum_Registered_PREVENTAD_2025-03-21.csv"),
  show_col_types = FALSE
)

# PET participant info (for dates and metadata)
pet_info <- read_csv(
  file.path(PREVENTAD_DATA_PATH,
            "PET_info_participants_Registered_PREVENTAD_2024-08-01.csv"),
  show_col_types = FALSE
)

# PROCESS PET DATA ------------------------------------------------------------

# Select key variables from NAV PET
# amyloid_index_SUVR = composite amyloid SUVR
# WhlCbl_centiloid = Centiloid scale (standardized)
pet_amyloid <- pet_nav %>%
  filter(PET_tracer == "NAV") %>%
  # Apply inclusion criteria
  filter(CONP_ID %in% criteria_met$CONP_ID) %>%
  select(CONP_ID, session, amyloid_index_SUVR, WhlCbl_centiloid) %>%
  rename(
    pet_session = session,
    amyloid_suvr = amyloid_index_SUVR,
    centiloid = WhlCbl_centiloid
  )

# Get PET visit dates from info file
pet_dates <- pet_info %>%
  filter(PET_tracer == "NAV") %>%
  select(CONP_ID, PET_session, Date_PET_visit, Age_PET_visit) %>%
  distinct() %>%
  rename(
    pet_session = PET_session,
    pet_date = Date_PET_visit,
    age_at_pet = Age_PET_visit
  )

# Merge PET data with dates
pet_amyloid <- pet_amyloid %>%
  left_join(pet_dates, by = c("CONP_ID", "pet_session"))

# DEFINE AMYLOID POSITIVITY ---------------------------------------------------

# Thresholds for amyloid positivity
CENTILOID_THRESHOLD_PRIMARY <- 22.32   # PREVENT-AD cohort-specific
CENTILOID_THRESHOLD_SECONDARY <- 24.4  # Standard threshold

pet_amyloid <- pet_amyloid %>%
  mutate(
    # Primary threshold (PREVENT-AD specific)
    amyloid_positive = as.integer(centiloid >= CENTILOID_THRESHOLD_PRIMARY),
    amyloid_status = factor(
      ifelse(centiloid >= CENTILOID_THRESHOLD_PRIMARY, "A+", "A-"),
      levels = c("A-", "A+")
    ),
    # Secondary threshold (for sensitivity analysis)
    amyloid_positive_secondary = as.integer(centiloid >= CENTILOID_THRESHOLD_SECONDARY),
    amyloid_status_secondary = factor(
      ifelse(centiloid >= CENTILOID_THRESHOLD_SECONDARY, "A+", "A-"),
      levels = c("A-", "A+")
    )
  )

# Summary
cat("\n=== AMYLOID PET SUMMARY ===\n")
cat("Total PET scans:", nrow(pet_amyloid), "\n")
cat("Unique participants:", n_distinct(pet_amyloid$CONP_ID), "\n")
cat("\nPrimary threshold (Centiloid >=", CENTILOID_THRESHOLD_PRIMARY, "):\n")
print(table(pet_amyloid$amyloid_status))
cat("\nSecondary threshold (Centiloid >=", CENTILOID_THRESHOLD_SECONDARY, "):\n")
print(table(pet_amyloid$amyloid_status_secondary))
cat("\nCentiloid distribution:\n")
print(summary(pet_amyloid$centiloid))

# SAVE PET DATA ---------------------------------------------------------------

saveRDS(pet_amyloid,
        file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_PET_amyloid.rds"))

cat("\nPET data saved to:",
    file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_PET_amyloid.rds"), "\n")
