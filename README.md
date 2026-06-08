# Using AD Genetic and Clinical Risk Burden Reports as Predictors of Blood-based Biomarker Levels and Trajectories

## Project Overview
This project examines how cumulative Alzheimer's disease (AD) risk burden—spanning genetic predisposition, family history, and modifiable lifestyle factors—relates to plasma pTau217 levels and the timing of tau positivity in cognitively normal older adults. The analysis uses data from the PREVENT-AD cohort. A manuscript draft and an AAIC 2026 abstract are in progress.

### Key Analyses
- **Risk Scores**: CogDRisk, LIBRA2, and AD genetic risk scores (GRS); composite risk burden score (0–4)
- **Biomarkers**: Aβ42/40 ratio, GFAP, NFL, p-tau181, p-tau217, p-tau231, p-tau217/Aβ42
- **Aim 1**: Cross-sectional and longitudinal associations between risk burden and biomarker levels
- **Aim 2**: PET-anchored positivity thresholds, SILA-derived individual trajectories, and Cox models predicting age of tau positivity onset

---

## Directory Structure

### `/data`
Contains raw, intermediate, and cleaned data files organized by processing stage.

#### `/data/raw`
Raw data files downloaded from PREVENT-AD repositories.

**Data Sources:**
- **PREVENT-AD**: Data accessed from Box (AndrewsLab/data/PREVENT-AD/)
- `GCST90027158_buildGRCh38.tsv.gz`: Bellenguez 2022 GWAS summary statistics
- `41588_2022_1024_MOESM4_ESM.xlsx`: Supplementary SNP data from Bellenguez et al.

**Note**: Raw data is stored on Box and accessed via CloudStorage. This directory contains softlinks to the original data.

#### `/data/intermediate`
Processed datasets organized by stage:
- `PREVENTAD_raw.rds`: Raw loaded PREVENT-AD data
- `PREVENTAD_biomarkers.rds`: Filtered longitudinal biomarker data
- `PREVENTAD_crsfactors.rds` / `PREVENTAD_crsfactors_imputed.rds`: CRS factor datasets (non-imputed and missForest-imputed)
- `PREVENTAD_cogd.rds`, `PREVENTAD_cogd_original.rds`, `PREVENTAD_cogd_scored.rds`, `PREVENTAD_cogd_original_scored.rds`: CogDRisk variables and scores (imputed and original)
- `PREVENTAD_libra2.rds`, `PREVENTAD_libra2_original.rds`, `PREVENTAD_libra2_scored.rds`, `PREVENTAD_libra2_scored_original.rds`: LIBRA2 variables and scores (imputed and original)
- `PREVENTAD_APOE.rds`: APOE genotype data
- `PREVENTAD_GWAS.rds`: GWAS data
- `PREVENTAD_GRS.rds`: Genetic risk scores
- `PREVENTAD_fhx.rds`: Family history data
- `PREVENTAD_MCI.rds`: MCI status data
- `PREVENTAD_risk_dat_baseline.rds`: Baseline risk data
- `PREVENTAD_singlevisi_scores.rds`: Risk scores for participants with a single biomarker visit
- `PREVENTAD_lmm_dose_interact.rds`: LMM dose-response interaction model objects
- `PREVENTAD_negative_ptau_slopes.rds`: Participants with negative observed pTau slopes (flagged for SILA sensitivity)
- `PREVENTAD_SILA_individual_fits.rds`: Individual SILA trajectory fits
- `Bellenguez_snps.rds`: Processed SNP list from Bellenguez GWAS
- `PREVENTAD_biomarkers_gwas_all.rds`: Biomarker data merged with GWAS IDs

#### `/data/cleaned`
Final analysis-ready datasets:
- `PREVENTAD_baseline.rds`: Baseline biomarkers + all predictors
- `PREVENTAD_lastvisit.rds`: Last-visit biomarkers + all predictors
- `PREVENTAD_longitudinal.rds`: All longitudinal biomarkers + predictors
- `PREVENTAD_predictors.rds`: Predictor variables only (age, sex, BMI, education, APOE, GRS, family history, CogDRisk, LIBRA2)

---

### `/workflow`

#### `/workflow/DATA_PIPELINE.md`
Documentation of the complete data processing pipeline from raw data to analysis-ready datasets, including design principles and key outputs at each step.

#### `/workflow/analysis/00_load_analysis_data.R`
Central script that combines all predictor sources with biomarker data to produce the four final analysis datasets in `/data/cleaned/`.

#### `/workflow/scripts`
Core scripts for data processing and functions.

**Configuration:**
- `config.R`: Central configuration file defining paths, packages, and constants

**Data Wrangling (run in order):**
- `data_wrangling/01_load_PREVENT-AD_data.R`: Load raw PREVENT-AD data from Box
- `data_wrangling/02_clean_PREVENT-AD.R`: Clean, filter, and impute baseline data; creates `biomarkers.rds`, `crsfactors.rds`, APOE, GWAS, and family history datasets
- `data_wrangling/03_recode_PREVENT-AD_CogD.R`: Recode and calculate CogDRisk scores
- `data_wrangling/04_recode_PREVENT-AD_libra2.R`: Recode and calculate LIBRA2 scores
- `data_wrangling/05_PREVENT-AD_PET.R`: Process PET imaging data for amyloid positivity analysis
- `data_wrangling/05_PREVENT-AD_singlevisi_scores.R`: Compute risk scores for participants with a single biomarker visit

**Functions:**
- `functions/CRS.R`: Risk score calculation functions (`calculate_cogdrisk`, `calculate_libra2`)
- `functions/CogDrisk_recode.R`: CogDRisk variable recoding
- `functions/LIBRA_recode.R`: LIBRA variable recoding
- `functions/functions_plot.R`: Custom plotting functions

**Other Scripts:**
- `AD_GRS.R`: AD genetic risk score calculations
- `GRS.R`: General GRS functions

#### `/workflow/analysis/eda`
Exploratory data analysis and cohort characterization:
- `PREVENT-AD_completecase.qmd/.html`: Complete case analysis
- `PREVENT-AD_fhx.qmd/.html`: Family history exploratory analysis
- `PREVENTAD_missing_data_viz.qmd/.html`: Missing data visualization
- `PREVENTAD_APOE_availability.R`: APOE genotype data availability checks
- `PREVENTAD_biomarker_outliers.R`: Biomarker outlier identification
- `PREVENTAD_cohort_filtering_analysis.R`: Cohort inclusion/exclusion accounting
- `PREVENTAD_final_cohort_summary.R`: Final cohort characterization
- `PREVENTAD_missing_data_analysis.R`: Missing data summary
- `verify_imputation.R`: Imputation quality verification
- `IMPUTATION_STRATEGY.md`: Documentation of the imputation approach and rationale
- `tables_PREVENT-AD.R` / `tables_PREVENT-AD.Rmd`: Table generation scripts

#### `/workflow/analysis/aim1`
**Aim 1: Risk Burden and Biomarker Level Associations**

*Primary analysis scripts:*
- `aim1_riskburden.qmd/.html`: Risk burden score calculation and descriptives
- `aim1_linearregression.qmd`: Multivariable linear regression (individual predictors) and cross-sectional composite score associations; Hochberg-corrected per-predictor family (12 tests per family)
- `aim1_linearregression_s.qmd`: Sensitivity version of linear regression analyses (excluding MCI progressors)
- `aim1_lme.qmd`: Linear mixed-effects models (LME) for longitudinal biomarker trajectories; composite × time and individual predictor × time interactions; Hochberg-corrected per-predictor family

*Exploratory / archived:*
- `aim1_crosssectional.qmd/.html`: Earlier cross-sectional analysis (superseded by aim1_linearregression.qmd)
- `aim1_longitudinal.qmd`: Earlier longitudinal analysis (superseded by aim1_lme.qmd)
- `aim1_longitudinal_riskburden_original.qmd`: Longitudinal analysis using original (non-imputed) risk scores
- `aim1_xsectional_stratified.qmd`: Stratified analyses (APOE, sex, GRS status)
- `aim1_associations.qmd`: Predictor association plots
- `aim1_crs_cutoffs.qmd`: Clinical cutoff analyses for risk scores
- `aim1_sensitivity_outliers.qmd/.html`: Sensitivity analyses excluding biomarker outliers

#### `/workflow/analysis/aim2`
**Aim 2: PET Anchoring, SILA Trajectories, and Time to Tau Positivity**
- `PREVENT-AD_PET.qmd/.html`: Predicts amyloid PET positivity from plasma biomarkers (ptau181, ptau217) using ROC analysis; derives Youden-index optimal positivity cutoffs (ptau181: 6.42 pg/mL; ptau217: 3.06 pg/mL) anchored at Centiloid ≥ 22.32
- `PREVENT-AD_PET_tau.qmd/.html`: Tau PET sensitivity analysis with alternative Centiloid threshold (≥ 18)
- `PREVENT-AD_SILA.qmd/.html`: SILA (Sampled Iterative Local Approximation) applied to longitudinal ptau181 and ptau217 to estimate individual trajectories relative to positivity thresholds
- `PREVENT-AD_SILA_tau.qmd`: SILA sensitivity analysis with alternative tau thresholds
- `SILA_demo.qmd`: SILA method demonstration and parameter exploration
- `aim2_cox_sila.qmd/.html`: Cox proportional hazards models (with left truncation for delayed study entry) testing whether risk burden components predict age at estimated tau positivity
- `aim2_cox_sila_bootstrap.qmd/.html`: Bootstrap-based uncertainty quantification for Cox model results
- `bootstrap.R`: Bootstrap helper functions

---

### `/results`

#### `/results/EDA`
- `/results/EDA/figures/`: Biomarker trajectories, missing data patterns, predictor distributions, imputation comparisons
- `/results/EDA/tables/`: Cohort summary tables, exclusion accounting, outlier counts, APOE distribution, missing data summaries

#### `/results/Aim1`
- `/results/Aim1/figures/`: Regression plots, trajectory plots, dose-response plots, diagnostic plots, sensitivity forest plots
- `/results/Aim1/tables/`: Model result tables, sensitivity comparison CSVs, longitudinal dose-response results, SILA individual fits and summary
- `/results/Aim1/stratified/`: Stratified results by APOE ε4 status and GRS

#### `/results/Aim2`
- `/results/Aim2/figures/`: SILA trajectory plots, combined biomarker plots
- `/results/Aim2/tables/`: Cox model results (`PREVENTAD_cox_sila_results.csv`, `PREVENTAD_cox_sila_results_primary.csv`)
- `bootstrap_cache.rds`: Cached bootstrap samples

---

### `/docs`
Documentation and metadata: project documentation, data dictionaries, analysis protocols.

### `/resources`
Softlinks to original data sources, literature references, and external documentation.

### `/archive`
Archived code no longer in active use but retained for reference, including earlier EDA scripts and a prior cross-sectional analysis notebook.

### `/epi217`
Course project materials (`epi217_final_project.qmd/.html`) using PREVENT-AD data; retained for reference.

---

## Data Processing Workflow

### Complete Pipeline

```
01_load_PREVENT-AD_data.R         # Load raw data from Box
    ↓
02_clean_PREVENT-AD.R             # Clean, filter, impute → crsfactors, biomarkers, APOE, GWAS, GRS, fhx
    ↓
    ├─ 03_recode_PREVENT-AD_CogD.R    # → cogd_scored
    └─ 04_recode_PREVENT-AD_libra2.R  # → libra_scored
         ↓
    05_PREVENT-AD_PET.R               # → PET positivity thresholds
    05_PREVENT-AD_singlevisi_scores.R # → single-visit risk scores
         ↓
00_load_analysis_data.R           # → baseline, lastvisit, longitudinal, predictors
    ↓
Aim 1 analysis scripts (aim1_*.qmd)
Aim 2 scripts: PREVENT-AD_PET.qmd → PREVENT-AD_SILA.qmd → aim2_cox_sila.qmd
```

See `workflow/DATA_PIPELINE.md` for full details on inputs, outputs, and design principles.

### Imputation Architecture
- **Centralized imputation**: missForest (random forest-based single imputation) runs once in `02_clean_PREVENT-AD.R`
- **Shared imputed data**: All downstream risk score scripts use the same `crsfactors_imputed` dataset
- **Parallel originals**: Non-imputed versions are retained for sensitivity analyses
- **OOB imputation error**: NRMSE = 0.027 (continuous variables, excellent); PFC = 0.210 (categorical variables, moderate); score-level SD ratios: CogDRisk = 1.013, LIBRA2 = 1.007

---

## Key Variables

### Risk Scores and Predictors
- **CogDRisk**: Dementia risk score (clinical and lifestyle factors)
- **LIBRA2**: Lifestyle for BRAin health index (12 risk/protective factors)
- **AD-GRS**: Alzheimer's disease polygenic risk score (Bellenguez 2022 GWAS weights)
- **Composite risk burden**: 0–4 score (APOE ε4 carrier, top-quartile GRS, ≥2 first-degree relatives with AD, PXS >1 SD above mean)
- **APOE**: ε4 carrier status
- **Family history**: ≥2 first-degree relatives with sporadic AD

### Biomarkers
**Primary analysis (Aim 1 — 6 biomarkers):**
- **Amyloid**: Aβ42/40 ratio
- **Neurodegeneration**: GFAP (glial fibrillary acidic protein), NfL (neurofilament light)
- **Tau**: p-tau181, p-tau217
- **Derived**: p-tau217/Aβ42 ratio

**Additional (Aim 2 / exploratory):**
- p-tau231

- **Positivity thresholds** (PET-anchored, Youden index, Centiloid ≥ 22.32): p-tau181 = 6.42 pg/mL; p-tau217 = 3.06 pg/mL

### Transformations
- Log transformations: GFAP, NfL, p-tau181, p-tau217, p-tau217/Aβ42 ratio
- Untransformed: Aβ42/40 ratio

---

## Usage Notes

### Running the Analysis
1. Set up data access to Box CloudStorage
2. Run `workflow/scripts/config.R` to configure paths and load packages
3. Run the data wrangling pipeline (scripts 01–05) in order
4. Run `workflow/analysis/00_load_analysis_data.R` to create final analysis datasets
5. Render Quarto documents in `workflow/analysis/aim1/` and `workflow/analysis/aim2/`

For Aim 2, the recommended rendering order is:
1. `PREVENT-AD_PET.qmd` (derives positivity cutoffs)
2. `PREVENT-AD_SILA.qmd` (fits individual trajectories)
3. `aim2_cox_sila.qmd` (Cox models using SILA output)
4. `aim2_cox_sila_bootstrap.qmd` (bootstrap uncertainty)

### Requirements
- R packages: tidyverse, ggplot2, gtsummary, lme4, lmerTest, ggpubr, missForest, here, pROC, silaR, survival, patchwork
- Access to PREVENT-AD data on Box

### Important Notes
- All analyses filter to baseline (at or before first biomarker measurement)
- missForest (random forest-based single imputation) is used for missing CRS factor data; OOB error: NRMSE = 0.027, PFC = 0.210; score-level SD ratios are 1.013 (CogDRisk) and 1.007 (LIBRA2), indicating acceptable score-level fidelity for both
- SILA models require ≥2 biomarker visits; participants with only one visit are handled separately via `05_PREVENT-AD_singlevisi_scores.R`
- Cox models use left truncation (delayed entry) to correct for differential study entry ages


---

## Contact Information
For questions or additional information, please contact the AndrewsLab research team.
