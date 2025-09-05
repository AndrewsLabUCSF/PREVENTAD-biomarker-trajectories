
# SETUP -------------------------------------------------------------------
# Load configuration file
source('workflow/scripts/config.R')

## Path to data
PREVENTAD_DATA_PATH <- 'data/PREVENT-AD'



# DEMOGRAPHICS TABLE ------------------------------------------------------
# Load demographics files
demographics_file <- file.path(PREVENTAD_DATA_PATH, "Demographics_Registered_PREVENTAD_2025-03-21.csv")
ad8_file <- file.path(PREVENTAD_DATA_PATH, "AD8_Registered_PREVENTAD_2025-03-21.csv")
questionnaire_file <- file.path(PREVENTAD_DATA_PATH, "SelfReport_Behavioral_Questionnaires_Registered_PREVENTAD_2024-08-01.csv")

demographics <- read.csv(demographics_file)
ad8 <- read.csv(ad8_file)
questionnaire <- read.csv(questionnaire_file)

ad8_baseline <- ad8 %>%
  filter(str_detect(Visit_label, "BL"))

missing_ad8 <- demographics %>%
  anti_join(ad8_baseline, by="CONP_ID")

table1_gtsummary <- demographics %>% 
  left_join(ad8_baseline, by="CONP_ID") %>%
  select(Sex, Education_years, Education_level) %>%
  tbl_summary() %>%
  as_gt()

