# Using AD Genetic and Clinical Risk Burden Reports as Predictors of Blood-based Biomarker Levels and Trajectories

## Project Overview
This project examines the associations between Alzheimer's disease (AD) genetic and clinical risk burden scores and blood-based biomarker levels and trajectories. The analysis uses data from the PREVENT-AD cohort.

### Key Analyses
- **Risk Scores**: CogDrisk, LIBRA2, and AD genetic risk scores (GRS)
- **Biomarkers**: Aβ42/40 ratio, GFAP, NFL, p-tau181, p-tau217, p-tau217/Aβ42
- **Study Designs**: Cross-sectional and longitudinal associations

---

## Directory Structure

### `/data`
Contains raw, intermediate, and cleaned data files organized by processing stage.

#### `/data/raw`
Raw data files downloaded from PREVENT-AD repositories. These files serve as the starting point for the data analysis pipeline.

**Data Sources:**
- **PREVENT-AD**: Data accessed from Box (AndrewsLab/data/PREVENT-AD/)

**Note**: Raw data is stored on Box and accessed via CloudStorage. This directory may contain softlinks to the original data.

#### `/data/intermediate`
Processed data organized by cohort:
- `/data/intermediate/PREVENTAD/`: PREVENT-AD cleaned datasets, risk score data (CogDRisk, LIBRA2), and imputed datasets

#### `/data/cleaned`
Final cleaned and filtered datasets ready for analysis:
- `/data/cleaned/`: Combined analysis-ready datasets
- `/data/cleaned/filtered/`: Subset datasets filtered by specific criteria

---

### `/workflow`
All analysis code organized by processing stage and analytical aim.

#### `/workflow/scripts`
Core scripts for data processing and functions:

**Configuration:**
- `config.R`: Central configuration file defining paths, packages, and constants

**Data Wrangling:**
- `/workflow/scripts/data_wrangling/PREVENT-AD/`: PREVENT-AD data loading, cleaning, and risk score calculation scripts
  - `01_load_PREVENTAD_data.R`: Load raw PREVENT-AD data
  - `02_clean_PREVENTAD_data.R`: Clean and filter baseline data
  - `03_recode_PREVENTAD_CogD.R`: Calculate CogDRisk scores
  - `04_recode_PREVENTAD_libra.R`: Calculate LIBRA scores

**Functions:**
- `/workflow/scripts/functions/`: Custom functions for risk score calculation, plotting, and analysis
  - `CRS.R`: Cognitive risk score functions (calculate_cogdrisk, calculate_libra2)
  - `CogDrisk_recode.R`: CogDRisk variable recoding functions
  - `LIBRA_recode.R`: LIBRA variable recoding functions
  - `functions_plot.R`: Custom plotting functions

**Other Scripts:**
- `AD_GRS.R`: AD genetic risk score calculations
- `GRS.R`: General genetic risk score functions

#### `/workflow/analysis`
Analysis notebooks organized by research aim:

**Exploratory Data Analysis (EDA):**
- `/workflow/analysis/eda/`: Exploratory analyses and complete case analyses

**Aim 1: Risk Burden and Biomarker Associations:**
- `/workflow/analysis/aim1/`: Cross-sectional and longitudinal analyses
  - `aim1_riskburden.qmd`: Risk burden score calculation and descriptives
  - `aim1_crosssectional.qmd`: Cross-sectional associations between risk scores and biomarkers
  - `aim1_longitudinal.qmd`: Longitudinal trajectories of biomarkers by risk score
  - `aim1_xsectional_stratified.qmd`: Stratified analyses by demographics
  - `aim1_associations.qmd`: Association analyses
  - `aim1_crs_cutoffs.qmd`: Clinical cutoff analyses for risk scores

---

### `/results`
Output files organized by analysis type:

#### `/results/EDA`
Exploratory data analysis outputs:
- `/results/EDA/figures/`: Descriptive plots, distributions, missing data patterns
- `/results/EDA/tables/`: Summary statistics tables

#### `/results/Aim1`
Aim 1 analysis outputs:
- `/results/Aim1/figures/`: Regression plots, diagnostic plots, biomarker trajectories
- `/results/Aim1/tables/`: Model results, association tables
- `/results/Aim1/stratified/`: Stratified analysis results by APOE, sex, etc.

---

### `/docs`
Documentation and metadata:
- Project documentation
- Data dictionaries
- Analysis protocols
- Manuscript drafts

---

### `/resources`
Additional resources and references:
- Softlinks to original data sources
- Literature references
- External documentation

---

### `/archive`
Archived code and analyses no longer in active use but retained for reference.

---

## Data Processing Workflow

### PREVENT-AD Pipeline
1. **Load raw data** (`01_load_PREVENTAD_data.R`)
2. **Clean and filter to baseline** (`02_clean_PREVENTAD_data.R`)
3. **Calculate CogDRisk scores** (`03_recode_PREVENTAD_CogD.R`)
4. **Calculate LIBRA scores** (`04_recode_PREVENTAD_libra.R`)
5. **Perform EDA and complete case analysis**
6. **Run cross-sectional and longitudinal analyses** (Aim 1)

### Imputation Architecture
- **Centralized imputation**: missForest runs once in the cleaning script 
- **Shared imputed data**: Both risk score scripts use the same `baseline_imputed` dataset
- **Benefits**: Ensures consistency, improves efficiency, maintains single source of truth

---

## Key Variables

### Risk Scores
- **CogDRisk**: Dementia risk score including clinical and lifestyle factors
- **LIBRA2**: LIfestyle for BRAin health index (12 risk/protective factors)
- **AD-GRS**: Alzheimer's disease genetic risk score

### Biomarkers
- **Amyloid**: Aβ42/40 ratio
- **Neurodegeneration**: GFAP (glial fibrillary acidic protein), NFL (neurofilament light)
- **Tau**: p-tau181, p-tau217, p-tau231
- **Derived**: p-tau217/Aβ42 ratio

### Transformations
- Log transformations: GFAP, NFL, p-tau181, p-tau217, p-tau217/Aβ42 ratio
- Square root: p-tau231
- Untransformed: Aβ42/40 ratio

---

## Usage Notes

### Running the Analysis
1. Set up data access to Box CloudStorage
2. Run `workflow/scripts/config.R` to configure paths and load packages
3. Follow cohort-specific pipeline (see Data Processing Workflow above)
4. Render Quarto documents in `workflow/analysis/` for reports

### Requirements
- R packages: tidyverse, ggplot2, gtsummary, lme4, lmerTest, ggpubr, missForest, here
- Access to PREVENT-AD data on Box

### Important Notes
- All analyses filter to baseline (at or before first biomarker measurement)
- BrANCH uses missForest imputation for missing baseline data
- Imputation quality metrics (SD ratios) inform choice of original vs imputed data for analyses
- CogDRisk imputation quality: SD ratio = 0.823 (acceptable)
- LIBRA2 imputation quality: SD ratio = 1.654 (concerning - use original scores)

---

## Contact Information
For questions or additional information, please contact the AndrewsLab research team.