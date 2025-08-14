# GWAS Tool Suite for Galaxy

This repository contains a suite of five tools designed to perform a complete **Genome-Wide Association Study (GWAS)** workflow within the **Galaxy** platform.  
The tools take raw, synthetic GWAS data as input and produce a final analysis, including visualization and investigation of significant genetic markers.

Each tool consists of:
- A **Python script** for backend logic
- An **XML wrapper** for the Galaxy user interface

---

## Workflow Overview

The five tools are designed to be run **in sequence**, where the output of one tool serves as the input for the next.

### 1. SNP & Sample QC
Initial data cleaning step that:
- Filters the raw genotype data to remove low-quality SNPs and samples.
- Uses user-defined thresholds for **missingness** and **Minor Allele Frequency (MAF)**.
- Ensures quality and reliability for downstream analyses.

### 2. Allele Frequency Calculator
After data cleaning, this tool:
- Calculates genetic statistics for each SNP.
- Counts the number of samples with each genotype (`0`, `1`, `2`).
- Computes the **MAF**, essential for the association test.

### 3. Association Test (Logistic Regression)
The core of the GWAS:
- Performs **logistic regression** for each SNP.
- Tests association between genotype and binary phenotype (case/control).
- Controls for covariates such as **age** and **sex**.
- Produces a **p-value** for each SNP.

### 4. Manhattan Plot Generator
Visualization tool that:
- Generates a **Manhattan plot** from the association test results.
- Plots statistical significance (`-log10(p-value)`) of every SNP across the genome.
- Allows easy identification of significant genetic "hits".

### 5. Windowed LD Calculator
Investigation tool for significant hits:
- Focuses on a genomic region around a **focal SNP**.
- Calculates **Linkage Disequilibrium (LD)** between SNPs in the window.
- Outputs a heatmap showing SNP correlations.

---

## Python Libraries Used

- **pandas**
- **numpy**
- **statsmodels**
- **matplotlib**
- **seaborn**
