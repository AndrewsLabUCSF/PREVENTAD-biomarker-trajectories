
# Load configuration file
source('workflow/scripts/config.R')

# Load datasets
all_files <- list.files(PREVENTAD_DATA_PATH, full.names=TRUE)

file_names <- all_files[str_detect(all_files, 
                                   "Auditory|BP_Pulse_Weight|CDR|Clinical_diagnosis|Demographics|Medical_history|Genetics|GWAS|Lab_Registered|MCI|Med_use|Plasma|Behavioral")]

# Read to list
PREVENTAD_raw <- c()
for (file in file_names) {
  data <- read.csv(file)
  PREVENTAD_raw[[file]] <- data
}

# Name list elements
names(PREVENTAD_raw) <- c("auditory", "bp_pulse_weight", "CDR_FU", "diagnosis", 
                          "demographics", "CDR_BL", "medical_history", "genetics", 
                          "GWAS", "lab", "MCI", "meduse", "plasma_4plex", "plasma_ptau217", 
                          "questionnaire")

# Save extracted data
saveRDS(PREVENTAD_raw, 
        file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_raw.rds"))