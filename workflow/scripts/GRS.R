
# SETUP -------------------------------------------------------------------
library(readxl)

# Load configuration file
source('workflow/scripts/config.R')

# Load datasets
gwas_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_GWAS.rds"))
bellenguez_results <- read_xlsx("data/raw/41588_2022_1024_MOESM4_ESM.xlsx", 
                                sheet="Supplementary Table 32", skip=2)

gwas <- gwas_raw

GRSweights <- data.frame(
  rsid = bellenguez_results$rsid,
  beta = bellenguez_results$Beta,
  effect_allele = bellenguez_results$`Risk allele`,
  pval = bellenguez_results$`P value`
)



# GRS CALCULATION ---------------------------------------------------------
# Function for GRS calculation
calculate_grs <- function(geno_matrix, weights_df) {
  # Parse rsid and allele from column names (format: rsid_allele)
  geno_cols   <- colnames(geno_matrix)
  geno_rsid   <- sub("_[^_]+$", "", geno_cols)
  geno_allele <- sub(".*_", "", geno_cols)
  
  # Match on rsid
  weights_idx <- match(geno_rsid, weights_df$rsid)
  matched     <- !is.na(weights_idx)
  
  geno_subset    <- geno_matrix[, matched, drop = FALSE]
  weights_subset <- weights_df[weights_idx[matched], ]
  
  # Flip dosage where PREVENT-AD allele != Bellenguez risk allele
  needs_flip <- geno_allele[matched] != weights_subset$effect_allele
  geno_subset[, needs_flip] <- 2L - geno_subset[, needs_flip]
  
  # Calculate GRS as matrix multiplication
  grs <- as.numeric(geno_subset %*% weights_subset$beta)
  return(grs)
}

# Convert to matrix 
id_cols <- c("CONP_ID", "CONP_CandID") 
geno_matrix <- as.matrix(gwas_raw[, !names(gwas_raw) %in% id_cols])
mode(geno_matrix) <- "integer"

# Calculate GRS
GRS <- calculate_grs(geno_matrix, GRSweights) 
gwas <- cbind(gwas_raw, GRS) %>%
  relocate(GRS, .after=CONP_CandID)

saveRDS(gwas, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_GRS.rds"))
