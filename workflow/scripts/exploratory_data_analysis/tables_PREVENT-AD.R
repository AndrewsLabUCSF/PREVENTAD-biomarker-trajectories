
# SETUP -------------------------------------------------------------------
# Load configuration file
source('workflow/scripts/config.R')
source('workflow/scripts/exploratory_data_analysis/load_PREVENT-AD_data.R')


# CLINICAL TABLE -----------------------------------------------------------

table1_gtsummary <- clinical_dataset_cogdrisk %>% 
  mutate(Education_level = case_when(Education_level == "High" ~ "High (>11 years)",
                                     Education_level == "Middle" ~ "Middle (8-11 years)",
                                     Education_level == "Low" ~ "Low (<8 years)",
                                     TRUE ~ NA) %>% factor(),
         Hypertension = case_when(Hypertension == 0 ~ "No",
                                  Hypertension == 1 ~ "Yes",
                                  NA ~ NA) %>% factor(),
         High_cholesterol = case_when(High_cholesterol == 0 ~ "<6.5 mmol/liter",
                                      High_cholesterol == 1 ~ ">6.5 mmol/liter",
                                      TRUE ~ NA) %>% factor(),
         Depression = case_when(Depression == 0 ~ "No",
                                Depression == 1 ~ "Yes",
                                TRUE ~ NA) %>% factor(),
         Atrial_fibrillation = case_when(Atrial_fibrillation == 0 ~ "No",
                                         Atrial_fibrillation == 1 ~ "Yes",
                                         TRUE ~ NA) %>% factor(),
         TBI = case_when(TBI == 0 ~ "No",
                         TBI == 1 ~ "Yes",
                         TRUE ~ NA) %>% factor()) %>%
  select(-CONP_ID) %>%
  tbl_summary(
    by = Sex,
    type = list(Hypertension ~ "categorical",
                High_cholesterol ~ "categorical",
                Depression ~ "categorical",
                Atrial_fibrillation ~ "categorical",
                TBI ~ "categorical"),
    label = list(Education_years = "Education (years)",
                 Education_level = "Education (level)",
                 BMI_category = "BMI (category)",
                 High_cholesterol = "High cholesterol",
                 Atrial_fibrillation = "Atrial fibrillation",
                 Diabetes_treatment = "Diabetes"),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})")
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Sex**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous; n (%) for categorical") %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_gt() %>%
  tab_row_group(label = md("*Medical History*"), rows = 13:34) %>%
  tab_row_group(label = md("**Demographics**"), rows = 1:12) %>%
  tab_style(
    style = list(cell_text(weight = "bold"), 
                 cell_fill(color = "#f8f9fa")),
    locations = cells_row_groups()
  )

# Save table 1
gtsave(table1_gtsummary, filename = file.path(DATA_PATHS$output$tables, "PREVENT-AD_table1.png"))



# CLINICAL TABLE ----------------------------------------------------------
# Load clinical files
bp_pulse_weight_file <- file.path(PREVENTAD_DATA_PATH, 
                                  "BP_Pulse_Weight_Registered_PREVENTAD_2025-03-21.csv")
lab_file <- file.path(PREVENTAD_DATA_PATH,
                      "Lab_Registered_PREVENTAD_2025-03-21.csv")

bp_pulse_weight <- read.csv(bp_pulse_weight_file)
lab <- read.csv(lab_file)


# Create tables
bp_pulse_weight_baseline <- bp_pulse_weight %>%
  filter(str_detect(Visit_label, "BL"))

lab_baseline <- lab %>%
  filter(str_detect(Visit_label, "BL"))

table2_gtsummary <- demographics %>% 
  left_join(bp_pulse_weight_baseline, by="CONP_ID") %>%
  left_join(lab_baseline, by="CONP_ID") %>%
  select(Sex, Systolic_blood_pressure, Diastolic_blood_pressure,
         total_cholesterol_value, Height, Weight) %>%
  tbl_summary(
    by = Sex,
    label = list(Sex = "Sex",
                 Education_years = "Education (years)",
                 Education_level = "Education (level)")
  ) %>%
  as_gt()
table1_gtsummary




# LIFESTYLE TABLE ---------------------------------------------------------
# Create tables
missing_questionnaires <- demographics %>% #33
  anti_join(questionnaire, by="CONP_ID")

missing_smoking <- questionnaire %>%
  anti_join(table3_gtsummary, by="CONP_ID")  #15

# smoking var
table3_gtsummary <- demographics %>%
  left_join(questionnaire, by="CONP_ID") %>%
  select(CONP_ID, closest_study_visit_label, smoking_present) %>%
  mutate(smoking_present_label = case_when(
    smoking_present == 0 ~ "I have never",
    smoking_present == 1 ~ "I don't smoke but I have smoked occasionally",
    smoking_present == 2 ~ "I don't smoke but I used to smoke every day",
    smoking_present == 3 ~ "I smoke but not every day",
    smoking_present == 4 ~ "I smoke every day",
    TRUE ~ NA
  )) %>%
  filter(!is.na(smoking_present)) %>%
  select(smoking_present_label) %>%
  tbl_summary() %>%
  as_gt()
table3_gtsummary


social_gtsummary <- demographics %>%
  left_join(questionnaire, by="CONP_ID") %>%
  select(CONP_ID, closest_study_visit_label, starts_with("social_")) %>%
  filter(if_all(starts_with("social_"), ~!is.na(.))) %>%
  group_by(CONP_ID) %>%
  slice(1) %>%
  ungroup() %>%
  select(-CONP_ID, -closest_study_visit_label) %>%
  tbl_summary() %>%
  as_gt()
social_gtsummary


exercise_gtsummary <- questionnaire %>%
  select(CONP_ID, closest_study_visit_label, starts_with("exer_curr_")) %>%
  group_by(CONP_ID) %>%
  slice(1)
  filter(complete.cases(-c(CONP_ID, closest_study_visit_label)))

