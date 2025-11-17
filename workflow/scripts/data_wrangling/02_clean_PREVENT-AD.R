
# SETUP -------------------------------------------------------------------

# Load files and libraries
source('workflow/scripts/config.R')

library(missForest)
library(rxnorm)
library(pharm)

clinical_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_clinical_raw.rds"))
lifestyle_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_lifestyle_raw.rds"))
fhx_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_fhx_raw.rds"))


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
  mutate(SU_rxcui = find_approx_rxcui(SU_medication),
         trt_type = case_when(SU_rxcui %in% antihypertensives$RxCui & !(SU_rxcui %in% supplements$RxCui) ~ "antihypertensive",
                              SU_rxcui %in% antidiabetics$RxCui & !(SU_rxcui %in% supplements$RxCui) ~ "antidiabetic",
                              SU_rxcui %in% antihyperlipidemic$RxCui & !(SU_rxcui %in% supplements$RxCui) ~ "antihyperlipidemic",
                              SU_rxcui %in% antidepressants$RxCui & !(SU_rxcui %in% supplements$RxCui) ~ "antidepressants",
                              TRUE ~ NA),
         trt_yn = if_else(is.na(trt_type), 0, 1)) %>%
  pivot_wider(names_from=trt_type, values_from=trt_yn, names_glue="medusage_{trt_type}") %>%
  select(CONP_ID, starts_with("medusage_anti")) %>%
  group_by(CONP_ID) %>%
  fill(medusage_antihyperlipidemic:medusage_antidepressants, .direction="updown") %>%
  distinct() %>%
  replace(is.na(.), 0)

# Save medication dataset
saveRDS(meduse, file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_meduse.rds"))


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
         ext_relatives_AD = if_else(other_family_members_AD > 0, 1, 0),
         
         # earliest age of onset
         min_family_onset = pmin(
           father_onset_age,
           mother_onset_age,
           sibling_onset_age_1,  
           sibling_onset_age_2,
           sibling_onset_age_3,
           na.rm = TRUE
         ),
         
         # single comprehensive fhx variable
         fhx_burden = case_when(
           # 3 = high
           FDR_AD >= 3 ~ 3,
           FDR_AD >= 2 & !is.na(min_family_onset) & min_family_onset < 65 ~ 3,
           !is.na(min_family_onset) & min_family_onset < 60 ~ 3,
           
           # 2 = moderate
           FDR_AD >= 2 ~ 2,
           !is.na(min_family_onset) & min_family_onset >= 60 & min_family_onset < 65 ~ 2,
           
           # 1 = low
           FDR_AD == 1 ~ 1
         ) %>% as.factor(),
         
         FDRAD_1ormore = if_else(FDR_AD == 1, 1, 2),
         FDR_ext_AD = case_when(
           var1 == 1 & ext_relatives_AD == 0 ~ 1,
           var1 == 1 & ext_relatives_AD == 1 ~ 2,
           var1 == 2 & ext_relatives_AD == 0 ~ 3,
           var1 == 2 & ext_relatives_AD == 1 ~ 4,
           TRUE ~ NA
         ) %>% as.factor()
  ) %>%
  relocate(FDR_AD, .after=CONP_ID) %>%
  relocate(other_family_members_AD, .after=FDR_AD) %>%
  relocate(other_maternal_family_members_AD, .after=other_family_members_AD) %>%
  relocate(other_paternal_family_members_AD, .after=other_maternal_family_members_AD) %>%
  relocate(ext_relatives_AD, .after=FDR_AD) %>%
  relocate(min_family_onset, .after=ext_relatives_AD) %>%
  relocate(fhx_burden, .after=min_family_onset)

# Save dataset
saveRDS(dat_fhx, file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_fhx_dat.rds"))


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
saveRDS(clinical_imp_raw, file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_clinical_imp_raw.rds"))
saveRDS(lifestyle_imp_raw, file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_lifestyle_imp_raw.rds"))
