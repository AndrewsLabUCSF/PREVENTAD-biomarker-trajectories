
# SETUP -------------------------------------------------------------------

# Load files
source('workflow/scripts/config.R')

clinical_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_clinical_raw.rds"))
lifestyle_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_lifestyle_raw.rds"))
data_cogd <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_cogd_dat.rds"))

colSums(is.na(clinical_raw))
colSums(is.na(lifestyle_raw))
65/