
# SETUP -------------------------------------------------------------------

# Load files and libraries
source('workflow/scripts/config.R')

library(missForest)
library(rxnorm)
library(pharm)

PREVENTAD_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_dat.rds"))

biomarkers_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_biomarkers.rds"))
clinical_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_clinical_raw.rds"))
fhx_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_fhx_raw.rds"))
genetics_raw <- PREVENTAD_dat$genetics
gwas_raw <- PREVENTAD_dat$GWAS
lifestyle_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_lifestyle_raw.rds"))
mci_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_MCI_raw.rds"))


# APOE --------------------------------------------------------------------
apoe_dat <- genetics_raw %>%
  #select(CONP_ID, CONP_CandID, APOE) %>%
  mutate(
    # e4 carrier 
    genotype = str_replace_all(APOE, " ", ""), 
    allele1 = as.numeric(str_extract(genotype, "^\\d")),  # First number
    allele2 = as.numeric(str_extract(genotype, "\\d$")),  # Last number
    e4_count = (allele1 == 4) + (allele2 == 4),
    e4_status = case_when(
      e4_count == 0 ~ "Non-carrier",
      e4_count == 1 ~ "Heterozygous",
      e4_count == 2 ~ "Homozygous"
    ),
    e4_status = factor(e4_status, levels=c("Non-carrier", "Heterozygous", "Homozygous")),
    
    apoe = case_when(
      APOE == "3 2" ~ "e2+",
      APOE == "3 3" ~ "e3/e3",
      APOE %in% c("4 2", "4 3", "4 4") ~ "e4+",
      TRUE ~ NA),
    apoe = factor(apoe, levels=c("e3/e3", "e2+", "e4+"))
    ) %>%
  select(CONP_ID, CONP_CandID, apoe, e4_status)
  
# Save APOE dataset
saveRDS(apoe_dat, file.path(DATA_CLEANED_PATH$cleaned, "PREVENTAD_apoe.rds"))  
  

# BIOMARKERS --------------------------------------------------------------
biomarkers <- biomarkers_raw %>%
  select(CONP_ID, CONP_CandID, ends_with("_label"), 
         Date_taken, Candidate_Age, all_of(BIOMARKER_VARS)) %>%
  rename(
    age = Candidate_Age,
    ab_ratio = AB_ratio_simoa_4plex,
    gfap = GFAP_simoa_4plex,
    nfl = NFL_simoa_4plex,
    ptau181 = ptau181_simoa_UGOT,
    ptau217 = ptau217_simoa_UGOT,
    ptau231 = ptau231_simoa_UGOT
  ) %>%
  group_by(CONP_ID) %>%
  mutate(
    ptau217_ab42_ratio = ptau217/AB42_simoa_4plex,
    age = age/12,
    baseline_date = first(Date_taken),
    years = interval(ym(baseline_date), ym(Date_taken)) / years(1),
    baseline_age = first(age)
  ) %>%
  select(-AB42_simoa_4plex)

# Save biomarker dataset
saveRDS(biomarkers, file.path(DATA_CLEANED_PATH$cleaned, "PREVENTAD_biomarkers.rds"))


# MEDICATION --------------------------------------------------------------
# only need to look for antihypertensives, antidiabetics, antihyperlipidemic, and antidepressants
antihypertensives <- getRxCuiViaMayTreat("hypertension")
antidiabetics <- getRxCuiViaMayTreat("diabetes")
antihyperlipidemic <- getRxCuiViaMayTreat("hypercholesterolemia")
antidepressants <- getRxCuiViaMayTreat("depressive")
supplements <- data.frame(
  RxCui = c("4419", "21406", "1895", "318224", "11416", "6574"),
  name = c("fish oil", "coqenzyme q10", "calcium", "flax seed oil", "zinc", "magnesium")
)

meduse <- clinical_raw %>%
  select(CONP_ID, SU_medication) %>%
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

# Save medication dataset
saveRDS(meduse, file.path(DATA_CLEANED_PATH$cleaned, "PREVENTAD_meduse.rds"))


# FAMILY HISTORY ----------------------------------------------------------
dat_fhx <- fhx_raw %>%
  rowwise() %>%
  mutate(FDR_AD = sum(father_dx_ad_dementia, 
                      mother_dx_ad_dementia, 
                      sibling_dx_ad_dementia_count, na.rm=TRUE),
         other_family_members_AD = if_else(is.na(other_family_members_AD), 
                                           0, 
                                           other_family_members_AD),
         other_maternal_family_members_AD = if_else(is.na(other_maternal_family_members_AD), 
                                                    0, 
                                                    other_maternal_family_members_AD),
         other_paternal_family_members_AD = if_else(is.na(other_paternal_family_members_AD), 
                                                    0, 
                                                    other_paternal_family_members_AD),
         # binary: extensive fhx AD == 1; no extensive fhx AD == 0
         ER_AD = if_else(other_family_members_AD > 0, 1, 0),
         
         FDRAD_1ormore = if_else(FDR_AD == 1, 1, 2),
         
         # final family history variable
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
saveRDS(dat_fhx, file.path(DATA_CLEANED_PATH$cleaned, "PREVENTAD_fhx_dat.rds"))


# MCI ---------------------------------------------------------------------
mci <- mci_raw %>%
  # recode MCI subtype to unimpaired, MCI, and dementia labels
  mutate(mci_status = case_when(RC_MCI == 0 ~ 0,
                                RC_MCI_subtype >= 2 & RC_MCI_subtype <= 5 ~ 1,
                                RC_MCI_subtype == 6 ~ 2),
         mci_status = factor(mci_status,
                             levels=c(0, 1, 2),
                             labels=c("CN", "MCI", "Dementia")))



# MISSING DATA IMPUTATION -------------------------------------------------
# Names of variables with any NA
clinical_vars_to_imp <- names(clinical_raw)[sapply(clinical_raw, anyNA)]
clinical_vars_to_imp <- clinical_vars_to_imp[clinical_vars_to_imp != "PRN_medication"]
lifestyle_vars_to_imp <- names(lifestyle_raw)[sapply(lifestyle_raw, anyNA)]

# Subset only the variables with any NA and recode to numeric or factor for MF
clinical_imp <- clinical_raw %>%
  select(all_of(clinical_vars_to_imp)) %>%
  mutate_at(c("head_injury_hospitalized", "head_injury_severe",
              "diagnosed_impairment"), as.factor)

lifestyle_imp <- lifestyle_raw %>%
  select(all_of(lifestyle_vars_to_imp)) %>%
  mutate(smoking_present = as.factor(smoking_present)) %>%
  mutate(exer_curr_act2_days = as.numeric(exer_curr_act2_days)) %>%
  mutate(across(ends_with("_intensity"), as.factor)) %>%
  mutate(across(starts_with("social_"), as.factor)) 

clinical_imp_mf <- missForest(clinical_imp, verbose=TRUE, maxiter=10, ntree=100)
clinical_imp_mf_dat <- cbind(clinical_imp_mf$ximp, CONP_ID=clinical_raw$CONP_ID) %>%
  relocate(CONP_ID) 
lifestyle_imp_mf <- missForest(lifestyle_imp, verbose=TRUE, maxiter=10, ntree=100)
lifestyle_imp_mf_dat <- cbind(lifestyle_imp_mf$ximp, CONP_ID=lifestyle_raw$CONP_ID) %>%
  relocate(CONP_ID) %>%
  mutate(across(starts_with("social_"), as.numeric))

clinical_imp_raw <- clinical_raw
clinical_imp_raw[names(clinical_imp_mf_dat)] <- clinical_imp_mf_dat
lifestyle_imp_raw <- lifestyle_raw
lifestyle_imp_raw[names(lifestyle_imp_mf_dat)] <- lifestyle_imp_mf_dat

# Save cleaned datasets
saveRDS(clinical_imp_raw, file.path(DATA_CLEANED_PATH$cleaned, "PREVENTAD_clinical_imp_raw.rds"))
saveRDS(lifestyle_imp_raw, file.path(DATA_CLEANED_PATH$cleaned, "PREVENTAD_lifestyle_imp_raw.rds"))


# FILTERING ---------------------------------------------------------------
# Filter out participants with no GWAS data
apoe_filtered <- apoe_dat %>%
  filter(CONP_ID %in% gwas_raw$CONP_ID)

biomarkers_filtered <- biomarkers %>%
  filter(CONP_ID %in% gwas_raw$CONP_ID) 

clinical_imp_filtered_raw <- clinical_imp_raw %>%
  filter(CONP_ID %in% gwas_raw$CONP_ID)

fhx_filtered <- dat_fhx %>%
  filter(CONP_ID %in% gwas_raw$CONP_ID)

lifestyle_imp_filtered_raw <- lifestyle_imp_raw %>%
  filter(CONP_ID %in% gwas_raw$CONP_ID)

mci_filtered <- mci %>%
  filter(CONP_ID %in% gwas_raw$CONP_ID)

meduse_filtered <- meduse %>%
  filter(CONP_ID %in% gwas_raw$CONP_ID)

# Save filtered datasets
saveRDS(apoe_filtered, file.path(DATA_CLEANED_PATH$filtered, 
                                 "PREVENTAD_apoe_filtered.rds"))
saveRDS(biomarkers_filtered, file.path(DATA_CLEANED_PATH$filtered, 
                                       "PREVENTAD_biomarkers_filtered.rds"))
saveRDS(clinical_imp_filtered_raw, file.path(DATA_CLEANED_PATH$filtered, 
                                             "PREVENTAD_clinical_imp_filtered_raw.rds"))
saveRDS(fhx_filtered, file.path(DATA_CLEANED_PATH$filtered, 
                                "PREVENTAD_fhx_filtered.rds"))
saveRDS(lifestyle_imp_filtered_raw, file.path(DATA_CLEANED_PATH$filtered, 
                                             "PREVENTAD_lifestyle_imp_filtered_raw.rds"))
saveRDS(mci_filtered, file.path(DATA_CLEANED_PATH$filtered, 
                                "PREVENTAD_MCI_filtered.rds"))
saveRDS(meduse_filtered, file.path(DATA_CLEANED_PATH$filtered, 
                                   "PREVENTAD_meduse_filtered.rds"))


