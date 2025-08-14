import argparse
import pandas as pd
import numpy as np
import matplotlib
# Use a non-interactive backend for running on servers without a display
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    """
    Main function to generate a Manhattan plot and a top hits file from GWAS results.
    """
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Generate a Manhattan plot from association test results.")
    
    parser.add_argument('--results', required=True, help="Path to the association results file from Tool 3.")
    parser.add_argument('--pval_thresh', type=float, default=5e-8, help="Genome-wide significance threshold (e.g., 5e-8).")
    parser.add_argument('--top_n', type=int, default=20, help="Number of top hits to report.")

    # Output files
    parser.add_argument('--out_png', required=True, help="Output path for the Manhattan plot PNG image.")
    parser.add_argument('--out_hits', required=True, help="Output path for the top hits TSV file.")
    
    args = parser.parse_args()

    # --- Data Loading and Preparation ---
    try:
        # The input file from Tool 3 should have 'snp_id' as the index column.
        df = pd.read_csv(args.results, sep='\t', index_col='snp_id')
    except FileNotFoundError as e:
        print(f"Error: Results file not found - {e}")
        return

    # Drop rows with missing p-values, chromosome, or position
    df.dropna(subset=['p_value', 'chrom', 'pos'], inplace=True)
    df = df[pd.to_numeric(df['p_value'], errors='coerce').notna()]
    
    # Handle cases where p-value is 0, which would cause log10 to fail
    if (df['p_value'] == 0).any():
        # Replace 0 with a very small number to avoid infinity issues
        min_p_val = df[df['p_value'] > 0]['p_value'].min()
        df.loc[df['p_value'] == 0, 'p_value'] = min_p_val * 0.1

    # Calculate -log10(p-value)
    df['log_p'] = -np.log10(df['p_value'].astype(float))

    # Prepare for plotting: sort by chromosome and position
    # Convert 'chrom' to a categorical type to ensure natural sorting (1, 2, ..., 10, ..., 22, X)
    df['chrom'] = df['chrom'].astype('category')
    df['chrom_num'] = df['chrom'].cat.codes
    df = df.sort_values(['chrom_num', 'pos'])

    # Calculate cumulative position for the x-axis
    df['cumulative_pos'] = 0
    last_pos = 0
    chrom_breaks = []
    for i, chrom_group in df.groupby('chrom_num'):
        df.loc[chrom_group.index, 'cumulative_pos'] = chrom_group['pos'] + last_pos
        new_max = df.loc[chrom_group.index, 'cumulative_pos'].max()
        chrom_breaks.append(new_max)
        last_pos = new_max

    # --- Plot Generation ---
    plt.figure(figsize=(16, 8))
    
    # Get alternating colors for chromosomes
    colors = ['#002366', '#708090'] # Dark blue and slate gray
    
    # Plot each chromosome
    for i, (chrom_num, group) in enumerate(df.groupby('chrom_num')):
        plt.scatter(group['cumulative_pos'], group['log_p'], color=colors[i % 2], s=10, alpha=0.8)

    # Add significance line
    plt.axhline(y=-np.log10(args.pval_thresh), color='red', linestyle='--', linewidth=1)
    
    # --- Aesthetics ---
    plt.title('Manhattan Plot of GWAS Results', fontsize=20)
    plt.xlabel('Genomic Position', fontsize=15)
    plt.ylabel('-log10(p-value)', fontsize=15)
    
    # Set chromosome labels on the x-axis
    chrom_centers = df.groupby('chrom_num')['cumulative_pos'].mean()
    plt.xticks(chrom_centers.values, df.groupby('chrom_num')['chrom'].first(), rotation=60, fontsize=9)
    plt.tick_params(axis='x', which='major', length=0) # Hide x-axis ticks
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()

    # --- Save the plot ---
    plt.savefig(args.out_png, dpi=300, format='png')
    print(f"Manhattan plot saved to {args.out_png}")

    # --- Generate Top Hits Report ---
    top_hits = df.sort_values('p_value').head(args.top_n)
    top_hits = top_hits[['chrom', 'pos', 'p_value', 'log_p', 'maf', 'beta']]
    top_hits.to_csv(args.out_hits, sep='\t')
    print(f"Top {args.top_n} hits saved to {args.out_hits}")

if __name__ == '__main__':
    main()
