This code accompanies the publication: 

## Improving Mercury Exposure Estimates: A Meta-analysis of Methylmercury-to-Total Mercury Ratios in Canadian Traditional Foods

# Abstract

**Background:** Methylmercury (MeHg) exposure poses health risks to Indigenous Peoples in Canada, who often rely on fish and wildlife as part of their traditional diet. Risk assessments are frequently overly conservative because most available monitoring data report total mercury (tHg) rather than MeHg, creating uncertainty—particularly for foods in which the proportion of tHg present as MeHg is lower or highly variable.

**Objective:** To conduct a meta-analysis of paired methylmercury (MeHg) and total mercury (tHg) data in traditional foods, develop a comprehensive database, and generate conversion factors to improve estimation of MeHg concentrations from total mercury data.

**Methods:** We conducted a critical review of paired MeHg and tHg data from 31 published studies and two unpublished datasets. Species- and tissue-specific MeHg:tHg ratios were estimated using random-effects meta-analysis. Conservative conversion factors were derived from upper confidence limits with sensitivity analyses. Conversion factors were evaluated using dietary exposure data from the Inuit Health Survey and the First Nations Food, Nutrition, and Environment Study. Variability in MeHg:tHg ratios was explored using meta-regression, covariate analyses, and within-species assessments.

**Results:** MeHg:tHg ratios were quantified for 65 species–tissue combinations, revealing strong and systematic tissue- and species-specific patterns. Conversion factors were derived for 38 species–tissues, with most values below 100%, particularly for liver and other internal organs, indicating that assuming total mercury equals MeHg frequently overestimates exposure. Applying these factors to dietary survey data substantially reduced overestimation of MeHg intake in most regions, while regions dominated by low-concentration foods showed limited benefit. Exploratory analyses indicated that MeHg:tHg ratios varied by tissue, species, and geographic location, with additional influences of age and size in some taxa.

**Significance:** This study provides a systematic database and validated conversion factors for MeHg in commonly consumed traditional foods. Findings improve hazard identification and risk assessment while underscoring the continued need for site-specific measurements and research to address remaining data gaps.

**Impact Statement:** This study advances mercury risk assessment for Indigenous Peoples by providing conversion factors to better estimate methylmercury (MeHg) concentrations in traditional foods when only total mercury data are available. The findings reduce uncertainty in exposure assessments, improve hazard characterization accuracy, and highlight key factors driving variability, thereby strengthening both public health protection and culturally relevant food safety evaluations.

**Keywords:** country food; exposure modeling; food safety; risk assessment; risk characterization

---

## MeHg : THg Conversion and Meta-Analysis Pipeline

## Overview

This repository contains a complete analytical pipeline for:

1. Converting mercury concentration data from dry weight (DW) to wet weight (WW)  
2. Estimating MeHg:THg ratios using meta-analysis  
3. Evaluating tissue-, species-, and group-level patterns  
4. Assessing covariate effects (age, sex, size)  
5. Analyzing spatial (latitudinal) patterns  
6. Testing robustness of species-grouped MeHg:THg conversion factors  

The workflow is designed for use with compiled mercury databases containing total mercury (THg), methylmercury (MeHg), moisture content, and associated biological metadata.

---

## MeHg : THg Conversion and Meta-Analysis Pipeline

## Overview

This repository contains a complete analytical pipeline for:

1. Converting mercury concentration data from dry weight (DW) to wet weight (WW)  
2. Estimating MeHg:THg ratios using meta-analysis  
3. Evaluating tissue-, species-, and group-level patterns  
4. Assessing covariate effects (age, sex, size)  
5. Analyzing spatial (latitudinal) patterns  
6. Testing robustness of species-grouped MeHg:THg conversion factors  

The workflow is designed for use with compiled mercury databases containing total mercury (THg), methylmercury (MeHg), moisture content, and associated biological metadata.

---

## Workflow and Script Order

The analytical pipeline is organized into sequential scripts. Scripts should be run in the order listed below unless otherwise noted.

---

### Step 01 — Convert Dry Weight to Wet Weight  
**Script:** `convert-to-ww.R`

**Purpose:**
- Converts mercury concentrations reported on a dry-weight basis to wet-weight values
- Reconstructs missing standard deviations (SD) or standard errors (SE) when possible
- Estimates missing moisture content using tissue- and species-specific means
- Produces a harmonized dataset suitable for downstream meta-analysis

**Key Processes:**
- SD ↔ SE conversion using reported sample sizes
- Moisture content imputation (species × tissue, then tissue-only)
- Dry-weight to wet-weight conversion for THg and MeHg

**Outputs:**
- `df_converted.csv`
- `converted_hg_summary.csv`

---

### Step 02 — MeHg:THg Meta-Analysis  
**Script:** `mehg-meta-analysis.R`

**Purpose:**
- Quantifies MeHg:THg ratios using random-effects meta-analysis
- Evaluates tissue- and species-specific patterns
- Produces publication-ready summary tables and figures

**Key Processes:**
- Ratio-of-means (ROM) meta-analysis
- Species- and tissue-level subgroup analyses
- Exploratory data analysis (boxplots, forest-style plots)
- Heterogeneity assessment and confidence interval estimation

**Outputs:**
- `rom_results.csv`
- `species_summary.csv`
- `species_tissue_meta.csv`
- Publication figures (JPEG/TIFF)

---

### Step 03 — Spatial (Latitudinal) Analysis  
**Script:** `mehg-ratio-spatial.R`

**Purpose:**
- Evaluates geographic variation in MeHg:THg ratios
- Tests latitude as a moderator of mercury speciation

**Key Processes:**
- Meta-regression with centered latitude (and tissue interactions where applicable)
- Species-specific spatial analyses (e.g., Polar Bear, Walleye, Beluga Whale, Ringed Seal)
- Regression and bubble plots

**Outputs:**
- Latitude-based regression figures
- Multi-panel publication figures (e.g., TIFF for manuscripts)

---

### Step 04 — Covariate Effects on MeHg:THg Ratios  
**Script:** `ratio-covariates.R`

**Purpose:**
- Examines biological covariates influencing MeHg:THg ratios
- Assesses whether age, sex, body size, or growth metrics modify ratios

**Key Processes:**
- Correlation analyses (e.g., ratio vs age/weight/length)
- Linear models and ANCOVA
- Meta-analysis / subgroup meta-analysis where applicable
- Meta-regression (e.g., using ROM objects)

**Outputs:**
- Covariate effect plots (boxplots, scatterplots with regressions)
- Model summaries and diagnostic figures

---

### Step 05 — Sensitivity & Robustness Testing  
**Script:** `sensitivity-analysis.R`  
**Dependency:** Requires outputs from `mehg-meta-analysis.R`

**Purpose:**
- Tests robustness and stability of grouped MeHg:THg conversion factors
- Quantifies uncertainty and inter-species variability within groups

**Key Processes:**
- Random-effects meta-analysis with heterogeneity metrics (I², τ²)
- Prediction interval estimation
- Inter-species variability assessment (within-group spread)
- Leave-one-species-out (LOSO) sensitivity analysis

**Outputs:**
- `species_in_groupings.csv`
- `group_meta_stats.csv`
- `species_spread.csv`
- `loso_results.csv`
- `loso_summary.csv`
