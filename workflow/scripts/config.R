
# LOAD PACKAGES -----------------------------------------------------------
pacman::p_load(tidyr, dplyr, ggplot2, gtsummary, gt, here, tidyverse, data.table)


# FUNCTIONS ---------------------------------------------------------------
FUNCTIONS_PATH <- here("workflow", "scripts", "functions")
COGD_RECODE_FN <- file.path(FUNCTIONS_PATH, "CogDrisk_recode.R")
LIBRA_RECODE_FN <- file.path(FUNCTIONS_PATH, "LIBRA_recode.R")
CRS_FN <- file.path(FUNCTIONS_PATH, "CRS.R")




# OUTPUT DATA PATHS -------------------------------------------------------
EDA_OUTPUT_PATH <- list(
  figures = here("results", "EDA", "figures"),
  tables = here("results", "EDA", "tables")
)

DATA_INTERMEDIATE_PATH <- here("data", "intermediate")

DATA_OUTPUT_PATHS <- list(
  output = list(
    figures = here("results", "figures"),
    tables = here("results", "tables")
  ),
  data = list(
    intermediate = here("data", "intermediate"),
    cleaned = here("data", "cleaned")
  )
)
