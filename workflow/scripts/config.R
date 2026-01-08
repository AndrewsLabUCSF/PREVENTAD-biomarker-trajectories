
# LOAD PACKAGES -----------------------------------------------------------
pacman::p_load(tidyr, dplyr, ggplot2, gtsummary, gt, here, tidyverse, data.table,
               lubridate, ggpubr, jtools, lme4, forcats, glue, ggsci, lmerTest)


# FUNCTIONS ---------------------------------------------------------------
FUNCTIONS_PATH <- here("workflow", "scripts", "functions")

COGD_RECODE_FN <- file.path(FUNCTIONS_PATH, "CogDrisk_recode.R")
LIBRA_RECODE_FN <- file.path(FUNCTIONS_PATH, "LIBRA_recode.R")
CRS_FN <- file.path(FUNCTIONS_PATH, "CRS.R")
PLOT_FN <- file.path(FUNCTIONS_PATH, "functions_plot.R")


# BIOMARKER NAMES ---------------------------------------------------------
BIOMARKER_VARS <- c("AB_ratio_simoa_4plex", "AB42_simoa_4plex", "GFAP_simoa_4plex", "NFL_simoa_4plex",
                    "ptau181_simoa_UGOT", "ptau217_simoa_UGOT", "ptau231_simoa_UGOT")
BIOMARKERS <- c("ab_ratio", "gfap", "nfl", "ptau181", "ptau217", "ptau231",
                "ptau217_ab42_ratio")
BIOMARKERS_TRANS <- c("ab_ratio", "gfap_log", "nfl_log", "ptau181_log", 
                      "ptau217_log", "ptau231_sqrt", "ptau217_ab42_ratio_log")


# INPUT DATA PATH ---------------------------------------------------------
RAW_DATA_PATH <- file.path("~/Library/CloudStorage/Box-Box/AndrewsLab/data/PREVENT-AD/data/")


# OUTPUT DATA PATHS -------------------------------------------------------
EDA_OUTPUT_PATH <- list(
  figures = here("results", "EDA", "figures"),
  tables = here("results", "EDA", "tables")
)

AIM1_OUTPUT_PATH <- list(
  figures = here("results", "Aim1", "figures"),
  tables = here("results", "Aim1", "tables"),
  stratified = here("results", "Aim1", "stratified")
)

DATA_INTERMEDIATE_PATH <- here("data", "intermediate")

DATA_CLEANED_PATH <- list(
  cleaned = here("data", "cleaned"),
  filtered = here("data", "cleaned", "filtered")
)

DATA_OUTPUT_PATHS <- list(
  output = list(
    figures = here("results", "figures"),
    tables = here("results", "tables")
  )
)
