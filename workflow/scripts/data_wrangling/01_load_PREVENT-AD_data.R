
# SETUP -------------------------------------------------------------------

# Load configuration file
source('workflow/scripts/config.R')

# Load datasets
files <- list.files('data/raw/PREVENT-AD', full.names=TRUE)

PREVENTAD_dat <- c()
for (file in files) {
  data <- read.csv(file)
  PREVENTAD_dat[[file]] <- data
}

names(PREVENTAD_dat) <- c("AD8", "auditory", "bp_pulse_weight", "diagnosis", 
                          "demographics", "medical_history", "genetics", "GWAS", 
                          "lab", "meduse", "plasma_4plex", "plasma_ptau217", 
                          "questionnaire")

# view1 <- PREVENTAD_dat$meduse %>%
#   group_by(CONP_ID) %>%
#   slice(1)

# VARIABLE SELECTION ----------------------------------------------------
# Clinical
clinical_raw <- PREVENTAD_dat$demographics %>%
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
               select(CONP_ID, past_depression, past_atrial_fibrillation, treatment_diabetes,
                      treatment_hypertension, treatment_hyperlipidemia)),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$meduse) %>%
              select(CONP_ID, SU_medication, PRN_medication) %>%
              group_by(CONP_ID) %>%
              slice(1),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$questionnaire %>%
               select(CONP_ID, head_injury_hospitalized, head_injury_severe) %>%
               filter_at(vars(head_injury_severe, head_injury_hospitalized), any_vars(!is.na(.)))),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$auditory %>%
               select(CONP_ID, diagnosed_impairment) %>%
               group_by(CONP_ID) %>%
               slice(1)),
            by="CONP_ID")


# Family history
PREVENTAD_dat$family_history <- PREVENTAD_dat$demographics %>%
  select(CONP_ID, father_dx_ad_dementia:other_paternal_family_members_AD)

fhx_raw <- PREVENTAD_dat$family_history


# Lifestyle table
lifestyle_raw <- PREVENTAD_dat$demographics %>%
  select(CONP_ID) %>%
  left_join((PREVENTAD_dat$questionnaire %>%
               select(CONP_ID, smoking_present) %>%
               filter(!is.na(smoking_present)) %>%
               group_by(CONP_ID) %>%
               slice_tail()),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$questionnaire %>%
               select(CONP_ID, gds_score) %>%
               group_by(CONP_ID) %>%
               slice(1)),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$questionnaire %>%
              select(CONP_ID, epoch_score_currently) %>%
              filter(!is.na(epoch_score_currently)) %>%
              group_by(CONP_ID) %>%
              slice(1)),
            by="CONP_ID") %>%
  left_join((PREVENTAD_dat$questionnaire) %>%
              select(CONP_ID, pittsburgh_total_score) %>%
              filter(!is.na(pittsburgh_total_score)) %>%
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


# Genetics
genetics <- PREVENTAD_dat$genetics %>%
  mutate(
    # Clean up any spacing issues
    apoe_genotype = str_replace_all(APOE, " ", ""),  # Remove spaces
    
    # Extract individual alleles
    allele1 = as.numeric(str_extract(apoe_genotype, "^\\d")),  # First number
    allele2 = as.numeric(str_extract(apoe_genotype, "\\d$")),  # Last number
    
    # Count alleles (primary variable for analysis)
    apoe_e2_count = (allele1 == 2) + (allele2 == 2),
    apoe_e3_count = (allele1 == 3) + (allele2 == 3),
    apoe_e4_count = (allele1 == 4) + (allele2 == 4),
    
    # e4 carrier status (binary)
    apoe_e4_carrier = if_else(apoe_e4_count > 0, 1, 0),
    
    # e4 carrier status (labeled factor for plots)
    apoe_e4_status = factor(apoe_e4_count,
                            levels = 0:2,
                            labels = c("Non-carrier", "One e4", "Two e4")),
    
    # APOE category
    apoe_category = case_when(
      apoe_genotype == "32" ~ "e2+",
      apoe_genotype == "33" ~ "e3/e3",
      apoe_genotype %in% c("42", "43", "44") ~ "e4+",
      TRUE ~ NA
      )
  ) %>%
  select(-c(allele1, allele2))

apoe <- genetics %>%
  select(CONP_ID, starts_with("apoe_"))


# Biomarker table
biomarkers <- PREVENTAD_dat$plasma_4plex %>%
  full_join(PREVENTAD_dat$plasma_ptau217, 
            by=c("CONP_ID", "CONP_CandID", "Study_visit_label", "Visit_label", 
                 "Date_taken", "Candidate_Age")) %>%
  mutate(Candidate_Age = Candidate_Age/12) %>%
  arrange(CONP_ID)


# Save intermediate datasets
saveRDS(clinical_raw, 
        file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_clinical_raw.rds"))
saveRDS(lifestyle_raw, 
        file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_lifestyle_raw.rds"))
saveRDS(fhx_raw,
         file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_fhx_raw.rds"))
saveRDS(genetics,
        file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_genetics.rds"))
saveRDS(PREVENTAD_dat$GWAS,
        file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_GWAS.rds"))
saveRDS(apoe,
        file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_apoe.rds"))
saveRDS(biomarkers,
        file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_biomarkers.rds"))
saveRDS(PREVENTAD_dat, 
        file.path(DATA_OUTPUT_PATHS$data$intermediate, "PREVENTAD_dat.rds"))
