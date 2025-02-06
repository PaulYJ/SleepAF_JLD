## Overview

This repository contains scripts and files for our work entitled “*Data-driven sleep phenotyping links atrial fibrillation to non-hypopnea-related sleep disturbances*”. The following sections provide an overview of the contents and organization.

## Contents

### Data preprocessing 
- **Data availability:** All PSG data are available via the National Sleep Research Resource (https://sleepdata.org/datasets/shhs). Full PSG and clinical/covariate data are available for all interested parties pending completion of a Data Access and Use Agreement and Institutional Review Board approval, as outlined on the National Sleep Research Resource website.
- **preprocessing.m**: Performs PSG and clinical data preprocessing tasks including data cleaning, normalization, and feature extraction.  
- **readXML.m**: Utility script for parsing XML files provided by SHHS website and extracting structured data for high-resolution recordings .  

### Dimension reduction and Model training

- **PCA.m**: Executes Principal Component Analysis (PCA) for dimensionality reduction.
- **UMAP.m**: Applies UMAP (Uniform Manifold Approximation and Projection) for nonlinear dimensionality reduction.  

### Model training ###

- **SVM.m**: Builds and evaluates Support Vector Machine (SVM) models.  

### Differential abundance testing and clustering ###

- **Milo.R**: Performs differential abundance analysis on PSG data using the Milo framework.  
- **Miloplot.R**: Generates visualizations for Milo analysis results (e.g., UMAP plots, bee-swarm plots).  

### Statistical analysis

- **Logistic Regression.m**: Implements logistic regression models for classification tasks.  
- **longitudinal_cohort_analysis.m**: Analyzes longitudinal cohort data to track trends over time.  
- **Chi_square.m**: Conducts Chi-square tests for statistical analysis.  
- **sleep_stage_hyponea_analysis.m**: Focuses on sleep stage classification and hypopnea event detection.  

### Intermediate Files
- All files generated during script execution (e.g., `.mat`, `.csv`, `.RData`) are intermediate outputs. Delete these to rerun pipelines from scratch.  

## How to Use

1. **Dependencies**:  
   - MATLAB (v2021a or later) with Statistics and Machine Learning Toolbox.  
   - R (v4.2+) with packages: `miloR`, `ggplot2`, `seurat`.  

2. **Execution Order**:  
   - Start with `preprocessing.m` to prepare data.  
   - Run analysis scripts (e.g., `Chi_square.m`, `SVM.m`) as needed.  
   - Use  `Miloplot.R` and code sections in each script for visualization.  

3. **Data Access**:  
   - Follow NSRR guidelines at [sleepdata.org/datasets/shhs](https://sleepdata.org/datasets/shhs).

---

*For detailed instructions, refer to comments within each script.*  