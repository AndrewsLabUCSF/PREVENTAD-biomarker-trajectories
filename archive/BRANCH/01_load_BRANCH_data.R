
# SETUP -------------------------------------------------------------------

# Load configuration file
source('workflow/scripts/config.R')

# Load BrANCH dataset
all_dat <- read.csv(file.path(RAW_BRANCH_DATA_PATH, "shea_dataset_2025-02-20.csv"))

id_cols <- c("PIDN", "UnQID", "age_at_DCDate")

suffixes <- names(all_dat) %>%
  str_extract("\\.[^.]+$") %>%  # Extract everything after last dot
  str_remove("^\\.") %>%         
  unique() %>%
  na.omit() %>%
  setdiff("shea")

# Split into separate dataframes
datasets <- map(suffixes, function(suffix) {
  # Select columns ending with this suffix
  cols <- names(all_dat)[str_detect(names(all_dat), paste0("\\.", suffix, "$"))]
  
  # Create dataframe with ID columns first, then data columns
  all_dat %>%
    select(all_of(id_cols), all_of(cols)) %>%
    rename_with(~str_remove(.x, paste0("\\.", suffix, "$")), 
                -all_of(id_cols))  # Don't rename ID columns
})

# Name list elements
names(datasets) <- suffixes



# VARIABLE SELECTION ------------------------------------------------------

# Select elements
table_names <- c("general", "demographics", "diet", "cogntive_activity_scale",
                 "physical", "clinical_labs", "physical_activity_scale","social_network_index", 
                 "genetics", "health_history", "quanterix", "diagnosis", "diagnosis_latest", 
                 "janssen")

datasets_small <- datasets[names(datasets) %in% table_names]

# Add dob.shea to demographics table
if("demographics" %in% names(datasets_small)) {
  datasets_small$demographics <- datasets_small$demographics %>%
    left_join(
      all_dat %>% select(PIDN, DCDate, dob = dob.shea),
      by = "PIDN",
      relationship = "many-to-many"
    )
}

# Save intermediate dataset
saveRDS(datasets_small, 
        file.path(DATA_INTERMEDIATE_PATH$BRANCH, "BRANCH_raw.rds"))
