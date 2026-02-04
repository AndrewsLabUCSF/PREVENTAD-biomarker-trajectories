
# SETUP -------------------------------------------------------------------

# Load files and libraries
source('workflow/scripts/config.R')

library(missForest)
library(rxnorm)
library(pharm)

PREVENTAD_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_raw.rds"))


# DATA FILTERING CRITERIA -------------------------------------------------
## Inclusion criteria:
##   - CONP_ID with GWAS
##   - CONP_ID with at least 2 biomarker visits

gwas <- PREVENTAD_raw$GWAS

biomarkers <- PREVENTAD_raw$plasma_4plex %>%
  full_join(PREVENTAD_raw$plasma_ptau217) %>%
  arrange(CONP_ID, Date_taken) %>%
  select(CONP_ID, Visit_label, Date_taken, Candidate_Age, all_of(BIOMARKER_VARS)) %>%
  rename(
    age = Candidate_Age,
    ab_ratio = AB_ratio_simoa_4plex,
    ab42 = AB42_simoa_4plex,
    gfap = GFAP_simoa_4plex,
    nfl = NFL_simoa_4plex,
    ptau181 = ptau181_simoa_UGOT,
    ptau217 = ptau217_simoa_UGOT
  ) %>%
  group_by(CONP_ID) %>%
  mutate(
    ptau217_ab42_ratio = ptau217/ab42,
    age = age/12,
    baseline_date = first(Date_taken),
    years = interval(ym(baseline_date), ym(Date_taken)) / years(1),
    baseline_age = first(age)
  ) %>%
  select(-ab42)

criteria_met <- biomarkers %>%
  group_by(CONP_ID) %>%
  summarise(
    n_visits = n(),
    has_2plus_visits = n() >= 2
  ) %>%
  filter(has_2plus_visits) %>%
  filter(CONP_ID %in% gwas$CONP_ID)

cat("Participants meeting criteria:", nrow(criteria_met), "\n")


# RISK COMPONENT SPECIFIC DATASETS -----------------------------------------
## Steps:
##   1. Select variables
##   2. Filter
##   3. Create simple calculated fields, rename variables, etc.
##   4. Impute missing data

## Biomarkers ----
biomarkers_filtered <- biomarkers %>%
  filter(CONP_ID %in% criteria_met$CONP_ID)

# Save dataset
saveRDS(biomarkers_filtered, file.path(DATA_INTERMEDIATE_PATH$base, "PREVENTAD_biomarkers.rds"))


## Genetics ----
# APOE and GWAS
apoe_dat <- PREVENTAD_raw$genetics %>%
  filter(CONP_ID %in% criteria_met$CONP_ID) %>%
  mutate(
    genotype = str_replace_all(APOE, " ", ""), 
    allele1 = as.numeric(str_extract(genotype, "^\\d")),  # First number
    allele2 = as.numeric(str_extract(genotype, "\\d$")),  # Last number
    apoe_e4_count = (allele1 == 4) + (allele2 == 4),
    apoe = case_when(
      APOE == "3 2" ~ "e2+",
      APOE == "3 3" ~ "e3/e3",
      APOE %in% c("4 2", "4 3", "4 4") ~ "e4+",
      TRUE ~ NA),
    apoe = factor(apoe, levels=c("e3/e3", "e2+", "e4+"))
  ) %>%
  select(CONP_ID, apoe, apoe_e4_count)

gwas_filtered <- gwas %>% filter(CONP_ID %in% criteria_met$CONP_ID)

# Save datasets
saveRDS(apoe_dat, file.path(DATA_INTERMEDIATE_PATH$base, "PREVENTAD_APOE.rds"))
saveRDS(gwas_filtered, file.path(DATA_INTERMEDIATE_PATH$base, "PREVENTAD_GWAS.rds"))


## CRS factors ----
crs_factors <- PREVENTAD_raw$demographics %>%
  filter(CONP_ID %in% criteria_met$CONP_ID) %>%
  select(CONP_ID, Sex, Education_years, Height) %>%
  # Change Sex to factor
  mutate(Sex = factor(Sex, levels=c("Female", "Male"))) %>%
  left_join(
    # bp_pulse_weight
    (PREVENTAD_raw$bp_pulse_weight %>% 
       group_by(CONP_ID) %>%
       slice(1) %>%
       select(CONP_ID, Age=Candidate_Age, Systolic_blood_pressure, Diastolic_blood_pressure) %>%
       mutate(Age = Age/12)
    ), 
    by="CONP_ID") %>%
  left_join( 
    (PREVENTAD_raw$bp_pulse_weight %>%  # Get first non-NA weight 
       group_by(CONP_ID) %>%
       filter(!is.na(Weight)) %>%
       slice(1) %>%
       select(CONP_ID, Weight)
    ),
    by="CONP_ID") %>%
  # Lab; get first non-NA lab measurements
  left_join(
    (PREVENTAD_raw$lab %>%
       filter(!is.na(total_cholesterol_value)) %>%
       group_by(CONP_ID) %>%
       slice(1) %>%
       select(CONP_ID, total_cholesterol_value)
    ),
    by="CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$lab %>%
       filter(!is.na(LDL_value)) %>%
       group_by(CONP_ID) %>%
       slice(1) %>%
       select(CONP_ID, LDL_value)
    ),
    by="CONP_ID") %>%
  left_join(
    (PREVENTAD_raw$lab %>%
       select(CONP_ID, hba1c_value) %>%
       filter(!is.na(hba1c_value)) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by="CONP_ID") %>%
  # Medical history
  left_join(
    (PREVENTAD_raw$medical_history %>%
       select(CONP_ID, past_depression, past_atrial_fibrillation, treatment_diabetes,
              treatment_hypertension, treatment_hyperlipidemia) 
    ),
    by="CONP_ID") %>%
  # Head injury variables
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, head_injury_hospitalized, head_injury_severe) %>%
       filter_at(vars(head_injury_severe, head_injury_hospitalized), 
                 any_vars(!is.na(.)))
    ),
    by="CONP_ID") %>%
  # Geriatric depression scale
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, gds_score) %>%
       filter(!is.na(gds_score)) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by="CONP_ID") %>%
  # Pittsburgh score for insomnia
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, pittsburgh_total_score) %>%
       filter(!is.na(pittsburgh_total_score)) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by="CONP_ID") %>%
  # Auditory
  left_join(
    (PREVENTAD_raw$auditory %>%
       select(CONP_ID, diagnosed_impairment, subjective_hearing_impairment) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by="CONP_ID") %>%
  # Smoking
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, smoking_present) %>%
       filter(!is.na(smoking_present)) %>%
       group_by(CONP_ID) %>%
       slice(1)
    ),
    by="CONP_ID") %>%
  # Epoch score (cognitive activity)
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, epoch_score_currently) %>%
       filter(!is.na(epoch_score_currently)) %>%
       group_by(CONP_ID) %>%
       slice(1)),
    by="CONP_ID") %>%
  # Physical activity
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, starts_with("exer_curr_")) %>%
       select(!matches("category")) %>%
       filter(!is.na(exer_curr_act1_intensity)) %>%
       group_by(CONP_ID) %>%
       slice(1) %>%
       # Convert character columns to numeric before recoding
       mutate(across(starts_with("exer_curr_") & where(is.character), as.numeric)) %>%
       # Recode exercise NAs to 0 for participants who have earlier activities
       rowwise() %>%
       mutate(
         # For act2: if act1 has data and act2 is NA, set act2 to 0
         across(starts_with("exer_curr_act2_"),
                ~if_else(!is.na(exer_curr_act1_intensity) & is.na(.), 0, .)),
         # For act3: if act1 has data and act3 is NA, set act3 to 0
         across(starts_with("exer_curr_act3_"),
                ~if_else(!is.na(exer_curr_act1_intensity) & is.na(.), 0, .)),
         # For act4: if act1 has data and act4 is NA, set act4 to 0
         across(starts_with("exer_curr_act4_"),
                ~if_else(!is.na(exer_curr_act1_intensity) & is.na(.), 0, .)),
         # For act5: if act1 has data and act5 is NA, set act5 to 0
         across(starts_with("exer_curr_act5_"),
                ~if_else(!is.na(exer_curr_act1_intensity) & is.na(.), 0, .))
       ) %>%
       ungroup()
    ),
    by="CONP_ID") %>%
  # Social activity
  left_join(
    (PREVENTAD_raw$questionnaire %>%
       select(CONP_ID, starts_with("social_life_")) %>%
       group_by(CONP_ID) %>%
       filter_at(vars(starts_with("social_life_")), any_vars(!is.na(.))) %>%
       slice(1)
    ),
    by="CONP_ID")


### Medication ----
# Only need to look for antihypertensives, antidiabetics, antihyperlipidemic, and antidepressants
antihypertensives <- getRxCuiViaMayTreat("hypertension")
antidiabetics <- getRxCuiViaMayTreat("diabetes")
antihyperlipidemic <- getRxCuiViaMayTreat("hypercholesterolemia")
antidepressants <- getRxCuiViaMayTreat("depressive")
supplements <- data.frame(
  RxCui = c("4419", "21406", "1895", "318224", "11416", "6574"),
  name = c("fish oil", "coqenzyme q10", "calcium", "flax seed oil", "zinc", "magnesium")
)

meduse <- PREVENTAD_raw$meduse %>%
  filter(CONP_ID %in% criteria_met$CONP_ID) %>%  
  select(CONP_ID, SU_medication) %>%
  group_by(CONP_ID) %>%
  slice(1) %>%
  separate_rows(SU_medication, sep=";") %>%
  mutate(SU_medication = if_else(SU_medication == "other SU medication(s)", NA, SU_medication)) %>%
  rowwise() %>%
  mutate(
    SU_rxcui = find_approx_rxcui(SU_medication),
    trt_type = case_when(
      SU_rxcui %in% antihypertensives$RxCui & !(SU_rxcui %in% supplements$RxCui) ~ "antihypertensive",
      SU_rxcui %in% antidiabetics$RxCui & !(SU_rxcui %in% supplements$RxCui) ~ "antidiabetic",
      SU_rxcui %in% antihyperlipidemic$RxCui & !(SU_rxcui %in% supplements$RxCui) ~ "antihyperlipidemic",
      SU_rxcui %in% antidepressants$RxCui & !(SU_rxcui %in% supplements$RxCui) ~ "antidepressants",
      TRUE ~ NA
    ),
    trt_yn = if_else(is.na(trt_type), 0, 1)) %>%
  pivot_wider(names_from=trt_type, values_from=trt_yn, names_glue="medusage_{trt_type}") %>%
  select(CONP_ID, starts_with("medusage_anti")) %>%
  group_by(CONP_ID) %>%
  fill(medusage_antihyperlipidemic:medusage_antidepressants, .direction="updown") %>%
  distinct() %>%
  replace(is.na(.), 0)

# Merge meduse data back with clinical table
crs_factors <- crs_factors %>% left_join(meduse, by="CONP_ID")

# Save CRS dataset
saveRDS(crs_factors, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_crsfactors.rds"))


## Family history ----
fhx <- PREVENTAD_raw$demographics %>%
  filter(CONP_ID %in% criteria_met$CONP_ID) %>%
  select(CONP_ID, father_dx_ad_dementia:other_paternal_family_members_AD) %>%
  rowwise() %>%
  mutate(
    # Number of first degree relatives with AD
    FDR_AD = sum(father_dx_ad_dementia, 
                 mother_dx_ad_dementia, 
                 sibling_dx_ad_dementia_count, na.rm=TRUE),
    # Recode NAs to 0 in "other" columns
    other_family_members_AD = if_else(is.na(other_family_members_AD), 
                                      0, 
                                      other_family_members_AD),
    other_maternal_family_members_AD = if_else(is.na(other_maternal_family_members_AD), 
                                               0, 
                                               other_maternal_family_members_AD),
    other_paternal_family_members_AD = if_else(is.na(other_paternal_family_members_AD), 
                                               0, 
                                               other_paternal_family_members_AD),
    # Binary y/n (1/0) extended relatives with AD
    ER_AD = if_else(other_family_members_AD > 0, 1, 0),
    # Binary 1 or more FDR with AD
    FDRAD_1ormore = if_else(FDR_AD == 1, 1, 2),
    # Final family history variable
    family_history = case_when(
      FDRAD_1ormore == 1 & ER_AD == 0 ~ 1,
      FDRAD_1ormore == 1 & ER_AD == 1 ~ 2,
      FDRAD_1ormore == 2 & ER_AD == 0 ~ 3,
      FDRAD_1ormore == 2 & ER_AD == 1 ~ 4,
      TRUE ~ NA
    ),
    family_history = factor(family_history,
                            levels=c(1, 2, 3, 4),
                            labels=c("1 FDR", "1 FDR + ext",
                                     "2+ FDR", "2+ FDR + ext"))
  ) %>%
  relocate(FDR_AD, .after=CONP_ID) %>%
  relocate(other_family_members_AD, .after=FDR_AD) %>%
  relocate(other_maternal_family_members_AD, .after=other_family_members_AD) %>%
  relocate(other_paternal_family_members_AD, .after=other_maternal_family_members_AD) %>%
  relocate(ER_AD, .after=FDR_AD) %>%
  relocate(family_history, .after=ER_AD)

# Save dataset
saveRDS(fhx, file.path(DATA_INTERMEDIATE_PATH$base, "PREVENTAD_fhx.rds"))



# MISSING CRS FACTORS IMPUTATION ------------------------------------------
# Prepare data for missForest
crs_factors_mf <- crs_factors %>%
  # Remove act2-5 exercise variables 
  select(-matches("exer_curr_act[2-5]")) %>%
  # Convert binary/categorical variables to factors
  mutate(smoking_present = as.factor(smoking_present)) %>%
  mutate(exer_curr_act1_intensity = as.factor(exer_curr_act1_intensity)) %>%
  # mutate(exer_curr_act1_weeks = as.factor(exer_curr_act1_weeks)) %>%
  mutate(social_life_frequency_phone_calls = as.factor(social_life_frequency_phone_calls)) %>%
  mutate(across(c(starts_with("past_"), starts_with("treatment_"), starts_with("head_"),
                  ends_with("_impairment"), starts_with("social_life_frequency_visit"), 
                  starts_with("medusage")), 
                as.factor)) %>%
  select(-CONP_ID)

# Set CONP_ID as rownames to keep linkage
rownames(crs_factors_mf) <- crs_factors$CONP_ID

# Run missForest imputation
crs_mf <- missForest(crs_factors_mf, verbose=TRUE, maxiter=10, ntree=100)
cat("Imputation complete. OOB error:", crs_mf$OOBerror, "\n")

# Extract imputed data
crs_factors_imputed <- crs_mf$ximp %>%
  rownames_to_column("CONP_ID") %>%
  # Convert factor variables back to numeric for downstream use
  mutate(
    across(c(starts_with("past_"), starts_with("treatment_"), starts_with("head_"),
             starts_with("medusage_")),
           ~as.numeric(as.character(.))),
    # Keep diagnosed_impairment and subjective_hearing_impairment as numeric
    across(ends_with("_impairment"), ~as.numeric(as.character(.))),
    # Convert social life frequency variables back to numeric
    across(starts_with("social_life_frequency"), ~as.numeric(as.character(.))),
    # Keep smoking_present and exercise variables as numeric
    smoking_present = as.numeric(as.character(smoking_present)),
    across(contains("_intensity"), ~as.numeric(as.character(.))),
    across(contains("_weeks"), ~as.numeric(as.character(.))),
    across(contains("_hours"), ~as.numeric(as.character(.)))
  )

# Get original act2-5 columns
act2_5_original <- crs_factors %>%
  select(CONP_ID, matches("exer_curr_act[2-5]"))

# Add back act2-5 exercise columns with original values where available, 0 where missing
crs_factors_imputed <- crs_factors_imputed %>%
  left_join(act2_5_original, by = "CONP_ID") %>%
  # Replace NAs with 0 for act2-5 columns
  mutate(across(matches("exer_curr_act[2-5]"), ~replace_na(., 0)))

# Save imputed dataset
saveRDS(crs_factors_imputed, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_crsfactors_imputed.rds"))



# MCI PROGRESSION ---------------------------------------------------------
mci <- PREVENTAD_raw$MCI %>%
  filter(CONP_ID %in% criteria_met$CONP_ID)

# Save MCI dataset
saveRDS(mci, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_MCI.rds"))
