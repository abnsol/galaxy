import argparse
import pandas as pd
import numpy as np

def main():
    """
    Main function to run the SNP and Sample QC pipeline.
    """
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Filter GWAS data based on SNP and sample QC metrics.")

    # Input files
    parser.add_argument('--genotypes', required=True, help="Path to the genotypes file (e.g., genotypes.tsv).")
    parser.add_argument('--snp_annot', required=True, help="Path to the SNP annotation file (e.g., snp_annotation.tsv).")
    parser.add_argument('--pheno', required=False, help="Optional: Path to the phenotypes/covariates file.")

    # QC thresholds
    parser.add_argument('--maf_thresh', type=float, default=0.01, help="Minimum Minor Allele Frequency (MAF) threshold.")
    parser.add_argument('--snp_miss', type=float, default=0.05, help="Maximum missingness threshold for SNPs.")
    parser.add_argument('--sample_miss', type=float, default=0.1, help="Maximum missingness threshold for samples.")
    parser.add_argument('--missing_val', type=str, default='NA', help="String representing missing values in the genotype file.")

    # Output files
    parser.add_argument('--out_geno', required=True, help="Output path for filtered genotypes.")
    parser.add_argument('--out_snp', required=True, help="Output path for filtered SNP annotations.")
    parser.add_argument('--out_pheno', required=False, help="Output path for filtered phenotypes/covariates.")
    parser.add_argument('--out_report', required=True, help="Output path for the QC summary report.")
    parser.add_argument('--out_samples', required=True, help="Output path for the list of samples that passed QC.")

    args = parser.parse_args()

    # --- Data Loading ---
    try:
        geno_df = pd.read_csv(args.genotypes, sep='\t', index_col=0)

        try:
            snp_annot_df = pd.read_csv(args.snp_annot, sep='\t', index_col='snp_id')
        except ValueError:
            snp_annot_df = pd.read_csv(args.snp_annot, sep='\t', index_col=0)

        pheno_df = None
        if args.pheno:
            pheno_df = pd.read_csv(args.pheno, sep='\t', index_col='sample_id')

    except FileNotFoundError as e:
        print(f"Error: Input file not found - {e}")
        return

    # --- Initial State ---
    initial_snp_count = geno_df.shape[1]
    initial_sample_count = geno_df.shape[0]

    # --- QC Step 1: Handle Missing Values ---
    geno_df.replace(args.missing_val, np.nan, inplace=True)
    geno_df = geno_df.astype(float)


    # --- QC Step 2: SNP Missingness Filter ---
    snp_missingness = geno_df.isnull().sum() / initial_sample_count
    snps_to_keep_miss = snp_missingness[snp_missingness <= args.snp_miss].index
    geno_df = geno_df[snps_to_keep_miss]
    
    snps_removed_miss = initial_snp_count - len(snps_to_keep_miss)


    # --- QC Step 3: MAF Filter ---
    snps_in_annot = snp_annot_df.index.intersection(geno_df.columns)
    maf_filtered_snps = snp_annot_df.loc[snps_in_annot]
    snps_to_keep_maf = maf_filtered_snps[maf_filtered_snps['maf'] >= args.maf_thresh].index
    geno_df = geno_df[snps_to_keep_maf]

    snps_removed_maf = len(snps_to_keep_miss) - len(snps_to_keep_maf)


    # --- QC Step 4: Sample Missingness Filter ---
    current_snp_count = geno_df.shape[1]
    sample_missingness = geno_df.isnull().sum(axis=1) / current_snp_count
    samples_to_keep = sample_missingness[sample_missingness <= args.sample_miss].index
    geno_df = geno_df.loc[samples_to_keep]

    samples_removed = initial_sample_count - len(samples_to_keep)


    # --- Final Filtering of Annotation and Phenotype Data ---
    final_snps = geno_df.columns
    final_samples = geno_df.index

    snp_annot_df = snp_annot_df.loc[final_snps]
    if pheno_df is not None:
        pheno_df = pheno_df.loc[final_samples]


    # --- Generate QC Report ---
    report_content = f"""
    --- SNP & Sample QC Report ---

    Initial Data:
    - {initial_snp_count} SNPs
    - {initial_sample_count} Samples

    Filtering Steps:
    1. SNP Missingness Filter (Threshold > {args.snp_miss}):
       - Removed {snps_removed_miss} SNPs.
       - {len(snps_to_keep_miss)} SNPs remaining.

    2. Minor Allele Frequency (MAF) Filter (Threshold < {args.maf_thresh}):
       - Removed {snps_removed_maf} SNPs.
       - {len(snps_to_keep_maf)} SNPs remaining.
    
    3. Sample Missingness Filter (Threshold > {args.sample_miss}):
       - Removed {samples_removed} samples.
       - {len(samples_to_keep)} samples remaining.

    Final Data:
    - {geno_df.shape[1]} SNPs passed all filters.
    - {geno_df.shape[0]} Samples passed all filters.
    """

    # --- Writing Output Files ---
    with open(args.out_report, 'w') as f:
        f.write(report_content)

    geno_df.fillna(-1, inplace=True)
    geno_df = geno_df.astype(int)
    geno_df.replace(-1, args.missing_val, inplace=True)

    geno_df.to_csv(args.out_geno, sep='\t', index_label='sample_id')
    snp_annot_df.to_csv(args.out_snp, sep='\t', index_label='snp_id')
    
    pd.DataFrame(final_samples, columns=['sample_id']).to_csv(args.out_samples, sep='\t', index=False)

    if args.out_pheno and pheno_df is not None:
        pheno_df.to_csv(args.out_pheno, sep='\t', index_label='sample_id')

    print("QC process completed successfully.")

if __name__ == '__main__':
    main()
