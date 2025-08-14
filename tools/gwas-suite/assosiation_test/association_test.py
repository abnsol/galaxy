import argparse
import pandas as pd
import numpy as np
import statsmodels.api as sm

def run_association_test(df, phenotype_col, genotype_col, covariate_cols):
    """
    Runs a single logistic regression test for one SNP.
    
    Args:
        df (pd.DataFrame): DataFrame containing phenotype, genotype, and covariates.
        phenotype_col (str): The name of the phenotype column.
        genotype_col (str): The name of the genotype column (the current SNP).
        covariate_cols (list): A list of strings for the covariate column names.
        
    Returns:
        A dictionary with the regression results (beta, SE, z, p_value) or None if the model fails.
    """
    try:
        # Drop rows with any missing values for the specific columns in this test
        model_df = df[[phenotype_col, genotype_col] + covariate_cols].dropna()

        # Ensure there's variation in the genotype and phenotype for the model to run
        if model_df[genotype_col].nunique() < 2 or model_df[phenotype_col].nunique() < 2:
            return None

        # Define the model variables
        Y = model_df[phenotype_col]
        X = model_df[[genotype_col] + covariate_cols]
        X = sm.add_constant(X, has_constant='add')

        # Fit the logistic regression model
        logit_model = sm.Logit(Y, X)
        result = logit_model.fit(disp=0)  # disp=0 suppresses convergence messages

        # Extract results for the genotype column
        beta = result.params[genotype_col]
        se = result.bse[genotype_col]
        z_score = result.tvalues[genotype_col]
        p_value = result.pvalues[genotype_col]
        
        return {'beta': beta, 'se': se, 'z': z_score, 'p_value': p_value}

    except Exception as e:
        return None

def main():
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Perform a SNP-by-SNP association test using logistic regression.")
    
    # Input files
    parser.add_argument('--genotypes', required=True, help="Path to the QC'd genotypes file.")
    parser.add_argument('--pheno_cov', required=True, help="Path to the phenotype and covariates file.")
    parser.add_argument('--snp_annot', required=False, help="Optional: SNP annotation file with MAF for filtering.")

    # Model parameters
    parser.add_argument('--phenotype_name', type=str, default='phenotype', help="Name of the phenotype column.")
    parser.add_argument('--covariates', type=str, default='', help="Comma-separated list of covariate column names.")
    parser.add_argument('--min_maf', type=float, default=0.01, help="Minimum MAF for a SNP to be tested.")
    parser.add_argument('--missing_val', type=str, default='NA', help="String for missing genotypes.")

    # Output file
    parser.add_argument('--output', required=True, help="Output path for the association results.")
    
    args = parser.parse_args()

    # --- Data Loading and Merging ---
    try:
        # Genotypes: Samples are rows (index), SNPs are columns
        geno_df = pd.read_csv(args.genotypes, sep='\t', index_col=0)
        
        # Phenotypes/Covariates: Samples are rows (index)
        pheno_cov_df = pd.read_csv(args.pheno_cov, sep='\t', index_col=0)

    except FileNotFoundError as e:
        print(f"Error: Input file not found - {e}")
        return
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    # Replace missing genotype placeholder with NaN for calculations
    geno_df.replace(args.missing_val, np.nan, inplace=True)

    merged_df = pheno_cov_df.join(geno_df, how='inner')
    
    # --- Covariate and SNP list setup ---
    covariate_list = [cov.strip() for cov in args.covariates.split(',') if cov.strip()]
    # SNPs are the columns of the original genotype dataframe
    snp_list = geno_df.columns.tolist()
    
    # --- MAF Filtering (Optional) ---
    snps_to_test = set(snp_list)
    snp_annot_df = None
    if args.snp_annot:
        try:
            snp_annot_df = pd.read_csv(args.snp_annot, sep='\t', index_col=0)
            snp_annot_df.index.name = 'snp_id'
            
            # Ensure 'maf' column exists before trying to filter
            if 'maf' in snp_annot_df.columns:
                snps_passing_maf = snp_annot_df[snp_annot_df['maf'] >= args.min_maf].index
                snps_to_test = set(snp_list).intersection(snps_passing_maf)
                print(f"Testing {len(snps_to_test)} SNPs after MAF filtering.")
            else:
                print("Warning: 'maf' column not found in annotation file. Skipping MAF filter.")

        except (FileNotFoundError, KeyError) as e:
            print(f"Warning: Could not use SNP annotation for MAF filtering: {e}")
            # Create a dummy dataframe to prevent downstream errors if file is bad
            snp_annot_df = pd.DataFrame(index=snp_list)
    
    if snp_annot_df is None:
        snp_annot_df = pd.DataFrame(index=snp_list)


    # --- Run Association Tests ---
    all_results = []
    for snp_id in snps_to_test:
        # The genotype for the current SNP is a column in our merged dataframe
        result = run_association_test(merged_df, args.phenotype_name, snp_id, covariate_list)
        
        if result:
            result['snp_id'] = snp_id
            all_results.append(result)

    # --- Format and Save Output ---
    if not all_results:
        print("Warning: No association tests could be completed successfully.")
        with open(args.output, 'w') as f:
            f.write("snp_id\tchrom\tpos\tbeta\tse\tz\tp_value\tmaf\n")
        return

    results_df = pd.DataFrame(all_results).set_index('snp_id')

    # Add annotation data (chrom, pos, maf) to the results
    cols_to_add = [col for col in ['chrom', 'pos', 'maf'] if col in snp_annot_df.columns]
    if cols_to_add:
        final_df = snp_annot_df.loc[results_df.index, cols_to_add].join(results_df)
    else:
        final_df = results_df # No annotation data to add

    # Reorder columns for clarity
    final_cols_order = ['chrom', 'pos', 'beta', 'se', 'z', 'p_value', 'maf']
    final_df = final_df.reindex(columns=final_cols_order)
    
    final_df.to_csv(args.output, sep='\t', na_rep='NA', index_label='snp_id')
    
    print(f"Association testing completed. Results for {len(final_df)} SNPs saved.")

if __name__ == '__main__':
    main()
