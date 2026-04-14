
# SETUP -------------------------------------------------------------------
library(readxl)
library(ggplot2)

# Load configuration file
source('workflow/scripts/config.R')

# Load datasets
gwas_raw <- readRDS(file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_GWAS.rds"))
bellenguez_results <- read_xlsx("data/raw/41588_2022_1024_MOESM4_ESM.xlsx", 
                                sheet="Supplementary Table 32", skip=2) %>%
  janitor::clean_names() %>%
    mutate(id = paste0(rsid,  "_", risk_allele))

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
  grs <- as.numeric(geno_subset %*% weights_subset$weight)
  return(grs)
}

# Convert to matrix 
id_cols <- c("CONP_ID", "CONP_CandID") 
geno_matrix <- as.matrix(gwas_raw[, !names(gwas_raw) %in% id_cols])
mode(geno_matrix) <- "integer"

# GRS QC ------------------------------------------------------------------
matched_bellenguez_rsids <- intersect(bellenguez_results$rsid, GRSweights$rsid)
unmatched_bellenguez_rsids <- setdiff(bellenguez_results$rsid, GRSweights$rsid)
unweighted_snp_cols <- GRSweights %>%
  filter(is.na(weight)) %>%
  pull(snp)

cat("=== GRS QC ===\n")
cat("GWAS SNP columns:", ncol(geno_matrix), "\n")
cat("Bellenguez variants:", nrow(bellenguez_results), "\n")
cat("Matched Bellenguez rsids:", length(matched_bellenguez_rsids), "\n")
cat("Unmatched Bellenguez rsids:", length(unmatched_bellenguez_rsids), "\n")
cat("SNP columns with missing weights:", length(unweighted_snp_cols), "\n")

if (length(unmatched_bellenguez_rsids) > 0) {
  cat("Unmatched Bellenguez rsids:\n")
  cat(paste(unmatched_bellenguez_rsids, collapse = ", "), "\n")
}

if (length(unweighted_snp_cols) > 0) {
  warning(
    "Some genotype columns do not have weights: ",
    paste(unweighted_snp_cols, collapse = ", ")
  )
}

# Calculate GRS
GRS <- calculate_grs(geno_matrix, GRSweights) 
gwas <- cbind(gwas_raw, GRS) %>%
  relocate(GRS, .after=CONP_CandID) %>% 
  as_tibble()

saveRDS(gwas, file.path(DATA_INTERMEDIATE_PATH, "PREVENTAD_GRS.rds"))

cat("\nGRS summary:\n")
print(summary(gwas$GRS))

grs_figure_dir <- here::here("results", "Aim2", "figures")
dir.create(grs_figure_dir, recursive = TRUE, showWarnings = FALSE)

grs_violin <- gwas %>%
  ggplot(aes(x = "PREVENT-AD", y = GRS)) +
  geom_violin(fill = "#67A9CF", color = "#2166AC", alpha = 0.8, linewidth = 0.5) +
  geom_boxplot(width = 0.12, outlier.shape = 21, outlier.size = 1.2, fill = "white") +
  labs(
    title = "PREVENT-AD Genetic Risk Score Distribution",
    x = NULL,
    y = "Genetic risk score"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(grs_figure_dir, "PREVENTAD_GRS_violin.png"),
  plot = grs_violin,
  width = 4.5,
  height = 4.5,
  dpi = 300
)
