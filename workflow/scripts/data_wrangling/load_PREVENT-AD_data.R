
# SETUP -------------------------------------------------------------------

# Load configuration file
source('workflow/scripts/config.R')

# Load datasets
files <- list.files('data/PREVENT-AD', full.names=TRUE)

PREVENTAD_dat <- c()
for (file in files) {
  data <- read.csv(file)
  PREVENTAD_dat[[file]] <- data
}

names(PREVENTAD_dat) <- c("AD8", "bp_pulse_weight", "diagnosis", "demographics",  
                          "medical_history", "genetics", "lab", "plasma_4plex", 
                          "plasma_ptau217", "questionnaire")


# COGDRISK --------------------------------------------------------------
## VARIABLE SELECTION ---------------------------------------------------

# Clinical
clinical_cogd_raw <- PREVENTAD_dat$demographics %>%
  select(CONP_ID, Sex, Education_years, Height) %>%
  left_join((PREVENTAD_dat$bp_pulse_weight %>% 
               group_by(CONP_ID) %>%
               slice(1) %>%
               select(CONP_ID, Candidate_Age, Systolic_blood_pressure, Diastolic_blood_pressure)), 
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$bp_pulse_weight %>%
               group_by(CONP_ID) %>%
               filter(!is.na(Weight)) %>%
               slice(1) %>%
               select(CONP_ID, Weight)),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$lab %>%
               filter(!is.na(total_cholesterol_value)) %>%
               group_by(CONP_ID) %>%
               slice(1) %>%
               select(CONP_ID, total_cholesterol_value)),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$lab %>%
               filter(!is.na(LDL_value)) %>%
               group_by(CONP_ID) %>%
               slice(1) %>%
               select(CONP_ID, LDL_value)),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$lab %>%
               select(CONP_ID, hba1c_value) %>%
               filter(!is.na(hba1c_value)) %>%
               group_by(CONP_ID) %>%
               slice(1)),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$medical_history %>%
               select(CONP_ID, past_depression, past_atrial_fibrillation, treatment_diabetes)),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$questionnaire %>%
               select(CONP_ID, head_injury_hospitalized, head_injury_severe) %>%
               filter_at(vars(head_injury_severe, head_injury_hospitalized), any_vars(!is.na(.)))),
            by="CONP_ID")


# Family history
PREVENTAD_dat$family_history <- PREVENTAD_dat$demographics %>%
  select(CONP_ID, father_dx_ad_dementia:other_paternal_family_members_AD)


# Lifestyle table
lifestyle_cogd_raw <- PREVENTAD_dat$demographics %>%
  select(CONP_ID) %>%
  left_join((PREVENTAD_dat$questionnaire %>%
               select(CONP_ID, smoking_present) %>%
               filter(!is.na(smoking_present)) %>%
               group_by(CONP_ID) %>%
               slice_tail()),
            by="CONP_ID") %>%
  left_join(PREVENTAD_dat$questionnaire %>%
              select(CONP_ID, epoch_score_currently) %>%
              filter(!is.na(epoch_score_currently)) %>%
              group_by(CONP_ID) %>%
              slice(1),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$questionnaire %>%
               select(CONP_ID, starts_with("exer_curr_")) %>%
               select(!matches("category")) %>%
               filter(!is.na(exer_curr_act1_intensity)) %>%
               group_by(CONP_ID) %>%
               slice(1)),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$questionnaire %>%
               select(CONP_ID, starts_with("social_life_")) %>%
               group_by(CONP_ID) %>%
               filter_at(vars(starts_with("social_life_")), any_vars(!is.na(.))) %>%
               slice(1)),
            by="CONP_ID")


# LIBRA -------------------------------------------------------------------

libra_dataset <- clinical_cogd_raw %>%
  select(CONP_ID, Systolic_blood_pressure, Diastolic_blood_pressure, total_cholesterol_value, 
         LDL_value, Height, Weight, hba1c_value, treatment_diabetes, past_depression) %>%
  left_join((lifestyle_cogd_raw %>%
               select(CONP_ID, smoking_present)),
            by="CONP_ID") %>%
  
# Recoding
  mutate(
    Hypertension = if_else(
      Systolic_blood_pressure >= 140 | Diastolic_blood_pressure >= 90, 1, 0),
    Hypercholesterolemia = if_else(
      total_cholesterol_value >= 5.2 & LDL_value >= 4.2, 1, 0),
    BMI = Weight/(Height/100)^2,
    BMI_category = case_when(BMI < 18.5 ~ "Underweight",
                             BMI > 18.5 & BMI < 25 ~ "Normal",
                             BMI >= 25 & BMI < 30 ~ "Overweight",
                             BMI >= 30 ~ "Obese",
                             TRUE ~ NA) %>% factor(),
    Diabetes = case_when(treatment_diabetes > 0 ~ 1,
                         hba1c_value >= 0.065 ~ 1,
                         TRUE ~ 0),
    Depression = past_depression %>% factor(),
    Smoking = if_else(smoking_present >= 3, 1, 0)
    ) %>%
  relocate(Hypertension, .after=Diastolic_blood_pressure) %>%
  relocate(Hypercholesterolemia, .after=LDL_value) %>%
  relocate(BMI, .after=Weight) %>%
  relocate(BMI_category, .after=BMI) %>%
  relocate(Diabetes, .after=BMI_category) %>%
  relocate(Depression, .after=past_depression) %>%
  relocate(Smoking, .after=smoking_present) %>%
  select(-c(Systolic_blood_pressure, Diastolic_blood_pressure, total_cholesterol_value, 
            LDL_value, Height, Weight, hba1c_value, treatment_diabetes,
            past_depression, smoking_present))

