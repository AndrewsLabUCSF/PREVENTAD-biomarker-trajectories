
# SETUP -------------------------------------------------------------------

# Load files and libraries
source('workflow/scripts/config.R')

library(missForest)
library(rxnorm)
library(pharm)

clinical_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_clinical_raw.rds"))
lifestyle_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_lifestyle_raw.rds"))


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