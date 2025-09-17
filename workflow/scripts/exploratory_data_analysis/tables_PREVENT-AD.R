
# SETUP -------------------------------------------------------------------
# Load configuration file
source('workflow/scripts/config.R')
source('workflow/scripts/exploratory_data_analysis/load_PREVENT-AD_data.R')


# CLINICAL TABLE -----------------------------------------------------------

table1_gtsummary <- clinical_dataset_cogdrisk %>% 
  # Changing labels
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
# table1_gtsummary

# Save table 1
gtsave(table1_gtsummary, filename = file.path(DATA_PATHS$output$tables, "PREVENT-AD_table1.png"))


# Split into demographics and medical history tables for slides
table1.1_gtsummary <- clinical_dataset_cogdrisk %>% 
  # Changing labels
  select(Age, Sex, Education_years, Education_level, BMI, BMI_category) %>%
  mutate(Education_level = case_when(Education_level == "High" ~ "High (>11 years)",
                                     Education_level == "Middle" ~ "Middle (8-11 years)",
                                     Education_level == "Low" ~ "Low (<8 years)",
                                     TRUE ~ NA) %>% factor()) %>%
  tbl_summary(
    by = Sex,
    label = list(Education_years = "Education (years)",
                 Education_level = "Education (level)",
                 BMI_category = "BMI (category)"),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})")
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Sex**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous; n (%) for categorical") %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_gt() %>%
  tab_row_group(label = md("**Demographics**"), rows = 1:12) %>%
  tab_style(
    style = list(cell_text(weight = "bold"), 
                 cell_fill(color = "#f8f9fa")),
    locations = cells_row_groups()
  )
#table1.1_gtsummary

table1.2_gtsummary <- clinical_dataset_cogdrisk %>% 
  select(Sex, Hypertension, High_cholesterol, Depression, Atrial_fibrillation, Diabetes_treatment, TBI) %>%
  mutate(Hypertension = case_when(Hypertension == 0 ~ "No",
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
  tbl_summary(
    by = Sex,
    type = list(Hypertension ~ "categorical",
                High_cholesterol ~ "categorical",
                Depression ~ "categorical",
                Atrial_fibrillation ~ "categorical",
                TBI ~ "categorical"),
    label = list(High_cholesterol = "High cholesterol",
                 Atrial_fibrillation = "Atrial fibrillation",
                 Diabetes_treatment = "Diabetes"),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})")
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Sex**") %>%
  modify_footnote(all_stat_cols() ~ "n (%) for categorical") %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_gt() %>%
  tab_row_group(label = md("*Medical History*"), rows = 1:22) %>%
  tab_style(
    style = list(cell_text(weight = "bold"), 
                 cell_fill(color = "#f8f9fa")),
    locations = cells_row_groups()
  )
table1.2_gtsummary

gtsave(table1.2_gtsummary, filename = file.path(DATA_PATHS$output$tables, "PREVENT-AD_table1.2.png"), zoom=1.2)

# LIFESTYLE TABLES --------------------------------------------------------

table2_gtsummary <- lifestyle_dataset_cogdrisk %>% 
  left_join((clinical_dataset %>% 
               select(CONP_ID, Sex)),
            by="CONP_ID") %>%
  select(-CONP_ID, -starts_with("social_life_")) %>%
  mutate(Smoking = case_when(Smoking == 2 ~ "Current",
                             Smoking == 1 ~ "Former",
                             Smoking == 0 ~ "Never",
                             TRUE ~ NA),
         Cognitive_engagement = case_when(Cognitive_engagement == 3 ~ "Highest",
                                          Cognitive_engagement == 2 ~ "Middle",
                                          Cognitive_engagement == 1 ~ "Lowest",
                                          TRUE ~ NA)) %>%
  tbl_summary(
    by = Sex,
    label = list(epoch_score_currently = "Epoch score (last year)",
                 Cognitive_engagement = "Cognitive engagement",
                 light_minutes_week = "Light",
                 moderate_minutes_week = "Moderate",
                 heavy_minutes_week = "Heavy"),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})")
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Sex**") %>%
  modify_footnote(all_stat_cols() ~ "Mean (SD) for continuous; n (%) for categorical") %>%
  modify_footnote_body(
    footnote = "Levels determined using tertiles of Epoch score",
    columns = label,
    rows = variable == "Cognitive_engagement" & row_type == "label"
  ) %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_gt() %>%
  tab_row_group(label = md("*Physical activity (mins per week)*"), rows = 13:15) %>%
  tab_style(
    style = list(cell_text(weight = "bold"), 
                 cell_fill(color = "#f8f9fa")),
    locations = cells_row_groups()
  )

# table2_gtsummary

# Save table 2
gtsave(table2_gtsummary, filename = file.path(DATA_PATHS$output$tables, "PREVENT-AD_table2.png"))


table3_gtsummary <- lifestyle_dataset_cogdrisk %>% 
  left_join((clinical_dataset %>% 
               select(CONP_ID, Sex)),
            by="CONP_ID") %>%
  select(-CONP_ID) %>%
  select(starts_with("social_life_"), Sex) %>%
  mutate(social_life_frequency_activities = 
           case_when(social_life_frequency_activities == 1 ~ "Less than 5 times a year",
                     social_life_frequency_activities == 2 ~ "5-10 times a year",
                     social_life_frequency_activities == 3 ~ "About once a month",
                     social_life_frequency_activities == 4 ~ "2-3 times a month",
                     social_life_frequency_activities == 5 ~ "About once a week",
                     social_life_frequency_activities == 6 ~ "Several days a week",
                     social_life_frequency_activities == 7 ~ "Every day",
                     TRUE ~ NA),
         social_life_frequency_activities = 
           factor(social_life_frequency_activities,
                  levels = c("Less than 5 times a year",
                             "5-10 times a year",
                             "About once a month",
                             "2-3 times a month",
                             "About once a week",
                             "Several days a week",
                             "Every day"),
                  ordered = TRUE),
         social_life_frequency_visitors = 
           case_when(social_life_frequency_visitors == 1 ~ "Not at all in past month",
                     social_life_frequency_visitors == 2 ~ "Once in past month",
                     social_life_frequency_visitors == 3 ~ "2-3 times in past month",
                     social_life_frequency_visitors == 4 ~ "About once a week",
                     social_life_frequency_visitors == 5 ~ "Several days a week",
                     social_life_frequency_visitors == 6 ~ "Every day",
                     TRUE ~ NA),
         social_life_frequency_visitors =
           factor(social_life_frequency_visitors,
                  levels = c("Not at all in past month",
                             "Once in past month",
                             "2-3 times in past month",
                             "About once a week",
                             "Several days a week",
                             "Every day"),
                  ordered = TRUE),
         social_life_frequency_visits =
           case_when(social_life_frequency_visits == 1 ~ "Not at all in past month",
                     social_life_frequency_visits == 2 ~ "Once in past month",
                     social_life_frequency_visits == 3 ~ "2-3 times in past month",
                     social_life_frequency_visits == 4 ~ "About once a week",
                     social_life_frequency_visits == 5 ~ "Several days a week",
                     social_life_frequency_visits == 6 ~ "Every day",
                     TRUE ~ NA),
         social_life_frequency_visits =
           factor(social_life_frequency_visits,
                  levels = c("Not at all in past month",
                             "Once in past month",
                             "2-3 times in past month",
                             "About once a week",
                             "Several days a week",
                             "Every day"),
                  ordered = TRUE),
         social_life_frequency_phone_calls =
           case_when(social_life_frequency_phone_calls == 1 ~ "Not at all in past month",
                     social_life_frequency_phone_calls == 2 ~ "Once in past month",
                     social_life_frequency_phone_calls == 3 ~ "2-3 times in past month",
                     social_life_frequency_phone_calls == 4 ~ "About once a week",
                     social_life_frequency_phone_calls == 5 ~ "Several days a week",
                     social_life_frequency_phone_calls == 6 ~ "Every day",
                     TRUE ~ NA),
         social_life_frequency_phone_calls =
           factor(social_life_frequency_phone_calls,
                  levels = c("Not at all in past month",
                             "Once in past month",
                             "2-3 times in past month",
                             "About once a week",
                             "Several days a week",
                             "Every day"),
                  ordered = TRUE)) %>%
  tbl_summary(
    by = Sex,
    label = list(social_life_frequency_activities = "Over a year's time, about how often do you get together with friends or relatives?",
                 social_life_frequency_visitors = "During the past month, about how often have you had friends over to your home?",
                 social_life_frequency_visits = "About how often have you visited with friends at their homes during the past month?",
                 social_life_frequency_phone_calls = "About how often were you on the telephone with close friends or relatives during the past month?"),
    statistic = list(all_continuous() ~ "{mean} ({sd})")
  ) %>%
  modify_header(label = "**Question**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Sex**") %>%
  modify_footnote(all_stat_cols() ~ "n (%) for categorical") %>%
  bold_labels() %>%
  italicize_levels() %>%
  as_gt() %>%
  tab_style(
    style = list(cell_text(weight = "bold"), 
                 cell_fill(color = "#f8f9fa")),
    locations = cells_row_groups()
  )
# table3_gtsummary

# Save table 3
gtsave(table3_gtsummary, filename = file.path(DATA_PATHS$output$tables, "PREVENT-AD_table3.png"))



