import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_r2(geno1, geno2):
    """
    Calculates the squared correlation coefficient (r^2) between two genotype vectors.
    """
    # Drop samples that have missing values in either SNP
    valid_indices = geno1.notna() & geno2.notna()
    if valid_indices.sum() < 2: # Need at least 2 samples to calculate correlation
        return np.nan
        
    geno1_clean = geno1[valid_indices]
    geno2_clean = geno2[valid_indices]
    
    # Pearson correlation coefficient
    corr = np.corrcoef(geno1_clean, geno2_clean)[0, 1]
    
    return corr ** 2

def main():
    """
    Main function to calculate and visualize Linkage Disequilibrium (LD) in a genomic window.
    """
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Calculate LD (r^2) for SNPs in a genomic window.")
    
    # Input files
    parser.add_argument('--genotypes', required=True, help="Path to the QC'd genotypes file.")
    parser.add_argument('--snp_annot', required=True, help="Path to the SNP annotation file.")
    
    # Window parameters
    parser.add_argument('--focal_snp', required=True, help="The ID of the SNP to center the window on.")
    parser.add_argument('--window_kb', type=int, default=250, help="Window size in kilobases (kb) on each side of the focal SNP.")
    parser.add_argument('--min_maf', type=float, default=0.01, help="Minimum MAF for SNPs to be included in the calculation.")
    parser.add_argument('--missing_val', type=str, default='NA', help="String for missing genotypes.")

    # Output files
    parser.add_argument('--out_matrix', required=True, help="Output path for the LD matrix (TSV).")
    parser.add_argument('--out_png', required=True, help="Output path for the LD heatmap (PNG).")
    
    args = parser.parse_args()

    # --- Data Loading ---
    try:
        geno_df = pd.read_csv(args.genotypes, sep='\t', index_col=0)
        snp_annot_df = pd.read_csv(args.snp_annot, sep='\t', index_col=0)
        snp_annot_df.index.name = 'snp_id'
    except FileNotFoundError as e:
        print(f"Error: Input file not found - {e}")
        return
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    # --- Find SNPs in Window ---
    if args.focal_snp not in snp_annot_df.index:
        print(f"Error: Focal SNP '{args.focal_snp}' not found in the annotation file.")
        return

    focal_info = snp_annot_df.loc[args.focal_snp]
    focal_chrom = focal_info['chrom']
    focal_pos = focal_info['pos']
    
    window_bp = args.window_kb * 1000
    window_start = max(0, focal_pos - window_bp)
    window_end = focal_pos + window_bp

    # Filter annotation for SNPs within the window on the same chromosome
    window_snps_df = snp_annot_df[
        (snp_annot_df['chrom'] == focal_chrom) &
        (snp_annot_df['pos'] >= window_start) &
        (snp_annot_df['pos'] <= window_end) &
        (snp_annot_df['maf'] >= args.min_maf)
    ]
    
    snps_in_window = window_snps_df.index.intersection(geno_df.columns)
    
    if len(snps_in_window) < 2:
        print("Fewer than 2 SNPs in the specified window. Cannot calculate LD.")
        # Create empty files to satisfy Galaxy
        open(args.out_matrix, 'w').close()
        open(args.out_png, 'w').close()
        return

    # --- LD Calculation ---
    # Subset genotype data to only include SNPs in our window
    geno_window_df = geno_df[snps_in_window].copy()
    geno_window_df.replace(args.missing_val, np.nan, inplace=True)
    
    # Create an empty DataFrame to store r^2 values
    ld_matrix = pd.DataFrame(index=snps_in_window, columns=snps_in_window, dtype=float)

    # Calculate r^2 for each pair of SNPs
    for i in range(len(snps_in_window)):
        for j in range(i, len(snps_in_window)):
            snp1_id = snps_in_window[i]
            snp2_id = snps_in_window[j]
            
            if i == j:
                r2 = 1.0
            else:
                r2 = calculate_r2(geno_window_df[snp1_id], geno_window_df[snp2_id])
            
            ld_matrix.loc[snp1_id, snp2_id] = r2
            ld_matrix.loc[snp2_id, snp1_id] = r2 # Matrix is symmetric

    # --- Generate Outputs ---
    ld_matrix.to_csv(args.out_matrix, sep='\t')
    print(f"LD matrix saved to {args.out_matrix}")

    # Create the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(ld_matrix.astype(float), cmap='Reds', vmin=0, vmax=1)
    plt.title(f"LD Heatmap for Chromosome {focal_chrom}\nWindow around {args.focal_snp}")
    plt.tight_layout()
    plt.savefig(args.out_png, dpi=300, format='png')
    print(f"LD heatmap saved to {args.out_png}")

if __name__ == '__main__':
    main()
