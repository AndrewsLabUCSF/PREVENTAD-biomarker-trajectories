
# SETUP -------------------------------------------------------------------

# Load configuration file
source('workflow/scripts/config.R')

# Load GWAS dataset
gwas_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_GWAS.rds"))
snps <- fread("data/raw/GCST90027158_buildGRCh38.tsv.gz")

saveRDS(snps, file.path(DATA_INTERMEDIATE_PATH, "Bellenguez_snps.rds"))
