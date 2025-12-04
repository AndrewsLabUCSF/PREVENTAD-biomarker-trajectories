
# SETUP -------------------------------------------------------------------
library(readxl)

# Load configuration file
source('workflow/scripts/config.R')

# Load datasets
gwas_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_dat.rds"))$GWAS
apoe <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_apoe.rds"))
bellenguez_results <- read_xlsx("data/raw/41588_2022_1024_MOESM4_ESM.xlsx", 
                                sheet="Supplementary Table 32", skip=2)

gwas <- gwas_raw

GRSweights <- data.frame(
  snp = paste0(bellenguez_results$rsid, "_", bellenguez_results$`Risk allele`), 
  beta = bellenguez_results$Beta,
  effect_allele = bellenguez_results$`Risk allele`,
  pval = bellenguez_results$`P value`
)



# GRS CALCULATION ---------------------------------------------------------
# Function for GRS calculation
calculate_grs <- function(geno_matrix, weights_df) {
  # Ensure SNPs are in same order
  common_snps <- intersect(colnames(geno_matrix), weights_df$snp)
  
  # Subset and reorder
  geno_subset <- geno_matrix[, common_snps, drop = FALSE]
  weights_subset <- weights_df[match(common_snps, weights_df$snp), ]
  
  # Calculate GRS as matrix multiplication
  grs <- as.numeric(geno_subset %*% weights_subset$beta)
  return(grs)
}

# Convert to matrix 
id_cols <- c("CONP_ID", "CONP_CandID") 
geno_matrix <- as.matrix(gwas_raw[, !names(gwas_raw) %in% id_cols])
mode(geno_imp_matrix) <- "integer"

# Calculate GRS
GRS <- calculate_grs(geno_matrix, GRSweights) 
gwas <- cbind(gwas_raw, GRS) %>%
  relocate(GRS, .after=CONP_CandID)

saveRDS(gwas, file.path(DATA_OUTPUT_PATHS$data$cleaned, "PREVENTAD_GRS.rds"))
