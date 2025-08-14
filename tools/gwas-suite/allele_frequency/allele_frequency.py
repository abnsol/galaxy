import argparse
import pandas as pd
import numpy as np

def calculate_maf(row):
    """
    Calculates the Minor Allele Frequency (MAF) for a single SNP (row).
    """
    count_het = row['count_1']
    count_hom_alt = row['count_2']
    n_non_missing = row['n_non_missing']

    if n_non_missing < 10:
        return np.nan
    
    # Calculate minor allele count and total alleles
    minor_allele_count = count_het + (2 * count_hom_alt)
    total_alleles = 2 * n_non_missing
    
    if total_alleles == 0:
        return np.nan
        
    maf = minor_allele_count / total_alleles
    return maf

def main():
    """
    Main function to calculate allele frequencies from genotype data.
    """
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Calculate allele frequencies from a genotype file.")
    
    # Input files
    parser.add_argument('--genotypes', required=True, help="Path to the filtered genotypes file.")
    parser.add_argument('--snp_annot', required=False, help="Optional: Path to the SNP annotation file to add chrom/pos.")
    parser.add_argument('--missing_val', type=str, default='NA', help="String representing missing values.")

    # Output file
    parser.add_argument('--output', required=True, help="Output path for the allele frequency results file.")
    
    args = parser.parse_args()

    # --- Data Loading ---
    try:
        geno_df = pd.read_csv(args.genotypes, sep='\t', index_col=0)
    except FileNotFoundError as e:
        print(f"Error: Genotype file not found - {e}")
        return

    # --- QC and Data Preparation ---
    geno_df.replace(args.missing_val, np.nan, inplace=True)
    
    # Transpose the dataframe so that SNPs are rows and samples are columns, which is easier for iteration.
    geno_df = geno_df.T

    # --- Allele Frequency Calculation ---
    results = []
    for snp_id, row in geno_df.iterrows():
        # Count occurrences of each genotype (0, 1, 2)
        counts = row.value_counts()
        
        # Get the number of non-missing samples for this SNP
        n_non_missing = row.count()
        
        # Create a dictionary for the current SNP's results
        snp_result = {
            'snp_id': snp_id,
            'count_0': counts.get(0, 0),
            'count_1': counts.get(1, 0),
            'count_2': counts.get(2, 0),
            'n_non_missing': n_non_missing
        }
        results.append(snp_result)

    # Convert the list of dictionaries to a DataFrame
    results_df = pd.DataFrame(results)
    results_df.set_index('snp_id', inplace=True)

    # Calculate MAF for each SNP using the dedicated function
    results_df['maf'] = results_df.apply(calculate_maf, axis=1)

    # --- Merge with Annotation Data (Optional) ---
    if args.snp_annot:
        try:
            snp_annot_df = pd.read_csv(args.snp_annot, sep='\t', index_col=0)
            
            snp_annot_df = snp_annot_df[['chrom', 'pos']]
            final_df = snp_annot_df.join(results_df, how='inner')
        except FileNotFoundError:
            print(f"Warning: SNP annotation file '{args.snp_annot}' not found. Output will not include chrom/pos.")
            final_df = results_df
        except KeyError:
            print("Error: The SNP annotation file must contain 'chrom' and 'pos' columns.")
            final_df = results_df
    else:
        final_df = results_df

    output_columns = ['chrom', 'pos', 'count_0', 'count_1', 'count_2', 'n_non_missing', 'maf']
    final_columns = [col for col in output_columns if col in final_df.columns]
    final_df = final_df[final_columns]

    # --- Writing Output File ---
    final_df.to_csv(args.output, sep='\t', na_rep='NA', index_label='snp_id') 
    
    print("Allele frequency calculation completed successfully.")

if __name__ == '__main__':
    main()
