# GWAS Workflow Run Log

**Date:** August 14, 2025  
**Analyst:** Abenezer  

This log documents the execution of the **5-tool GWAS workflow** on the provided synthetic dataset.

---

## Step 1: Data Generation
- **Script:** `generate_gwas_data.py`
- **Purpose:** Create the initial dataset in the `gwas_data/` directory.
- **Files Created:**
  - `genotypes.tsv`
  - `phenotypes_covariates.tsv`
  - `snp_annotation.tsv`

---

## Step 2: Tool Execution in Galaxy
The tools were executed **in order**, using default parameters unless otherwise specified.

### 1. SNP & Sample QC
- **Inputs:** `genotypes.tsv`, `snp_annotation.tsv`, `phenotypes_covariates.tsv` (from Step 1)
- **Parameters:** Default  
  - MAF > 0.01  
  - SNP Missingness < 0.05  
  - Sample Missingness < 0.1  
- **Outputs:**
  - Filtered genotypes  
  - Filtered annotations  
  - Filtered phenotypes  
  - QC report  
- **Result:** The synthetic data passed all filters, so no SNPs or samples were removed.

---

### 2. Allele Frequency Calculator
- **Inputs:** Filtered genotypes, filtered annotations
- **Parameters:** Default
- **Output:** Table of allele frequencies and genotype counts for each SNP

---

### 3. Association Test (Logistic Regression)
- **Inputs:** Filtered genotypes, filtered phenotypes, allele frequency table (as SNP annotation)
- **Parameters:**
  - Covariates: `age`, `sex`
  - Other: Default
- **Output:** Association results table with p-values for each SNP

---

### 4. Manhattan Plot Generator
- **Inputs:** Association results table
- **Parameters:**  
  - Significance threshold: `5e-8`  
  - Top N hits: `20`  
- **Outputs:**
  - PNG Manhattan plot
  - TSV file of top 20 hits

---

### 5. Windowed LD Calculator
- **Inputs:** Filtered genotypes, allele frequency table (as SNP annotation)
- **Parameters:**
  - Focal SNP: `rs100999` (selected from top hits)
  - Window size: `250 kb`
- **Outputs:**
  - PNG LD heatmap
  - TSV LD matrix

---

## Notes
- The workflow completed **successfully** from start to finish.  
- The synthetic dataset produced a **clear**, but **not genome-wide significant**, signal in the association test.  
- The Manhattan plot correctly visualized the results.  
- The LD calculator correctly identified a **region of low LD** around the focal SNP.
