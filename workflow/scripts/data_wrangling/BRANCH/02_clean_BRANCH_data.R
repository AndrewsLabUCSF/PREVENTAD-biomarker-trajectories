
# SETUP -------------------------------------------------------------------

# Load files and libraries
source('workflow/scripts/config.R')

library(missForest)

BRANCH_dat <- readRDS(file.path(DATA_INTERMEDIATE_PATH$BRANCH, "BRANCH_raw.rds"))

# Dataframe list for cleaned data
BRANCH_cleaned_dat <- c()


# BIOMARKER ---------------------------------------------------------------
ptau217_df <- BRANCH_dat$janssen %>%
  group_by(PIDN) %>%
  filter(DCDate != "") %>%
  filter(sum(!is.na(ptau217_janssen)) >= 2) %>%
  ungroup() %>%
  select(PIDN, UnQID, DCDate, gfap_janssen, nfl_janssen, ptau217_janssen)

BRANCH_cleaned_dat$janssen <- ptau217_df

quanterix_df <- BRANCH_dat$quanterix %>%
  group_by(PIDN) %>%
  filter(DCDate != "") %>%
  filter(
    sum(!is.na(ab40_branch)) >= 2 | 
      sum(!is.na(ab42_branch)) >= 2 | 
      sum(!is.na(gfap_branch)) >= 2 |
      sum(!is.na(nfl_branch)) >= 2 |
      sum(!is.na(ptau181_branch)) >= 2 
  ) %>%
  filter(
    if_any(c(ab40_branch, ab42_branch, gfap_branch, 
             nfl_branch, ptau181_branch), 
           ~!is.na(.))
  ) %>%
  ungroup() %>%
  select(PIDN, UnQID, DCDate, 
         ab40 = ab40_branch, ab42 = ab42_branch, 
         gfap = gfap_branch,
         nfl = nfl_branch, 
         ptau181 = ptau181_branch) %>%
  mutate(ab_ratio = ab42/ab40)

BRANCH_cleaned_dat$quanterix <- quanterix_df


# APOE --------------------------------------------------------------------

df <- BRANCH_dat$genetics %>%
  filter(PIDN %in% ptau217_df$PIDN | PIDN %in% quanterix_df$PIDN) %>%
  group_by(PIDN) %>%
  slice(1) %>%
  select(PIDN, UnQID, ApoE) %>%
  filter(str_detect(ApoE, "E")) %>%
  mutate(
    apoe = case_when(
      grepl("E4", ApoE) ~ "e4+",
      ApoE == "E3/E3" ~ "e3/e3",
      grepl("E2", ApoE) ~ "e2+",
      TRUE ~ NA),
    apoe = factor(apoe, levels=c("e3/e3", "e2+", "e4+"))
  ) %>%
  select(-ApoE)

BRANCH_cleaned_dat$APOE <- df

# Filter out individuals in biomarker tables with no APOE data
BRANCH_cleaned_dat$quanterix <- BRANCH_cleaned_dat$quanterix %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN)

BRANCH_cleaned_dat$janssen <- BRANCH_cleaned_dat$janssen %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN)


# FIRST BIOMARKER DATE ----------------------------------------------------

# Calculate first biomarker date for baseline filtering
first_biomarker <- bind_rows(
  BRANCH_cleaned_dat$janssen %>% select(PIDN, DCDate),
  BRANCH_cleaned_dat$quanterix %>% select(PIDN, DCDate)
) %>%
  mutate(DCDate = as.Date(DCDate)) %>%
  group_by(PIDN) %>%
  summarise(first_biomarker_date = min(DCDate, na.rm = TRUE), .groups = "drop")


# DEMOGRAPHICS ------------------------------------------------------------

df <- BRANCH_dat$demographics %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN) %>%
  group_by(PIDN) %>%
  slice(1) %>%
  mutate(gender = if_else(gender == 2, "Female", "Male")) %>%
  select(PIDN, UnQID, baselineDate = DCDate, dob, age_at_DCDate, gender, educ)

BRANCH_cleaned_dat$demographics <- df


# DIAGNOSIS ---------------------------------------------------------------

diagnosis_df <- BRANCH_dat$diagnosis %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN) %>%
  group_by(PIDN) %>%
  filter(res_dx_a != "")

diagnosis_latest_df <- BRANCH_dat$diagnosis_latest %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN) %>%
  select(-age_at_DCDate) %>%
  group_by(PIDN) %>%
  filter(clin_syn_best_est != "") %>%
  slice(1)

BRANCH_cleaned_dat$diagnosis <- diagnosis_df
BRANCH_cleaned_dat$diagnosis_latest <- diagnosis_latest_df


# HEALTH HISTORY ----------------------------------------------------------

df <- BRANCH_dat$health_history %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN) %>%
  filter(DCDate != "") %>%
  mutate(DCDate = as.Date(DCDate)) %>%
  left_join(first_biomarker, by = "PIDN") %>%
  filter(DCDate <= first_biomarker_date) %>%
  arrange(PIDN, DCDate) %>%
  group_by(PIDN) %>%
  slice(1) %>%
  select(PIDN, UnQID, DCDate, age_at_DCDate, diabetes, dep2yrs, quitsmok, smokyrs,
         alcfreq, cvangio, cvbypass, cvhatt, cbstroke, tbi, cvafib, insomn) %>%
  ungroup()

BRANCH_cleaned_dat$health_history <- df  


# PHYSICAL ----------------------------------------------------------------

df <- BRANCH_dat$physical %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN) %>%
  filter(DCDate != "") %>%
  arrange(PIDN, DCDate) %>%
  group_by(PIDN) %>%
  summarise(
    UnQID = first(UnQID),
    DCDate_physical = first(DCDate),
    height = first(height[!is.na(height)]),
    weight = first(weight[!is.na(weight)]),
    bpsys = first(bpsys[!is.na(bpsys)]),
    bpdias = first(bpdias[!is.na(bpdias)]),
    .groups = "drop"
  )

BRANCH_cleaned_dat$physical <- df  



# CLINICAL LABS -----------------------------------------------------------

df <- BRANCH_dat$clinical_labs %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN) %>%
  filter(DCDate != "") %>%
  mutate(DCDate = as.Date(DCDate)) %>%
  left_join(first_biomarker, by = "PIDN") %>%
  filter(DCDate <= first_biomarker_date) %>%
  filter(!is.na(total_cholesterol_mg_dl)) %>%
  select(PIDN, UnQID, DCDate, first_biomarker_date, total_cholesterol_mg_dl) %>%
  arrange(PIDN, DCDate) %>%
  group_by(PIDN) %>%
  slice(1) 

BRANCH_cleaned_dat$clinical_labs <- df  


# PHYSICAL ACTIVITY SCALE -------------------------------------------------

df <- BRANCH_dat$physical_activity_scale %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN) %>%
  filter(DCDate != "") %>%
  mutate(DCDate = as.Date(DCDate)) %>%
  left_join(first_biomarker, by = "PIDN") %>%
  filter(DCDate <= first_biomarker_date) %>%
  arrange(PIDN, DCDate) %>%
  group_by(PIDN) %>%
  slice(1) %>%
  select(PIDN, UnQID, DCDate, pase_total = pase_pase_total) %>%
  ungroup()

BRANCH_cleaned_dat$physical_activity_scale <- df  


# SOCIAL NETWORK INDEX ----------------------------------------------------

df <- BRANCH_dat$social_network_index %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN) %>%
  filter(DCDate != "") %>%
  mutate(DCDate = as.Date(DCDate)) %>%
  left_join(first_biomarker, by = "PIDN") %>%
  filter(DCDate <= first_biomarker_date) %>%
  arrange(PIDN, DCDate) %>%
  group_by(PIDN) %>%
  slice(1) %>%
  select(PIDN, UnQID, DCDate_sni = DCDate, 
         sni_network_size = SNI_NumberPeople, 
         sni_high_contact = SNI_HighContact) %>%
  ungroup()

BRANCH_cleaned_dat$social_network_index <- df  


# COGNITIVE ACTIVITY ------------------------------------------------------

df <- BRANCH_dat$cogntive_activity_scale %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN) %>%
  filter(DCDate != "") %>%
  mutate(DCDate = as.Date(DCDate)) %>%
  left_join(first_biomarker, by = "PIDN") %>%
  filter(DCDate <= first_biomarker_date) %>%
  arrange(PIDN, DCDate) %>%
  group_by(PIDN) %>%
  slice(1) %>%
  select(PIDN, UnQID, DCDate, cas_tot_pts) %>%
  ungroup()

BRANCH_cleaned_dat$cogntive_activity_scale <- df  


# DIET --------------------------------------------------------------------

df <- BRANCH_dat$diet %>%
  filter(PIDN %in% BRANCH_cleaned_dat$APOE$PIDN) %>%
  filter(DCDate != "") %>%
  mutate(DCDate = as.Date(DCDate)) %>%
  left_join(first_biomarker, by = "PIDN") %>%
  filter(DCDate <= first_biomarker_date) %>%
  arrange(PIDN, DCDate) %>%
  group_by(PIDN) %>%
  slice(1) %>%
  select(PIDN, UnQID, DCDate, dietq_mindscore) %>%
  ungroup()

BRANCH_cleaned_dat$diet <- df  


# Save cleaned dataset
saveRDS(BRANCH_cleaned_dat, file.path(DATA_INTERMEDIATE_PATH$BRANCH, "BRANCH_cleaned.rds"))
