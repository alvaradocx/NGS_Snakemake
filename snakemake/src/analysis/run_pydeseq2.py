#!/usr/bin/env python3
"""
run_pydeseq2.py
PyDESeq2 differential expression analysis script for multiple conditions
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from itertools import combinations

# PyDESeq2 imports
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# anndata/scanpy 
import scanpy as sc
import anndata as ad

# Plotting
from sklearn.decomposition import PCA
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

def parse_arguments():
    parser = argparse.ArgumentParser(description='PyDESeq2 differential expression analysis')
    parser.add_argument('--counts', required=True, help='Count matrix file')
    parser.add_argument('--samples', required=True, help='Sample metadata info file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--padj', type=float, default=5e-8, help='Adjusted p-value threshold')
    parser.add_argument('--lfc', type=float, default=1.0, help='Log2 fold change threshold')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--comparisons', help='Specific comparisons to run (format: cond1_vs_cond2,cond3_vs_cond4). If not specified, all pairwise comparisons will be performed.')
    
    return parser.parse_args()

def load_data(counts_file, meta_file):
    """Load count matrix and sample information"""
    print("Loading count matrix and sample information...")
    
    # Load counts and metadata info
    counts_df = pd.read_csv(counts_file, sep='\t', index_col=0)
    metadata = pd.read_csv(meta_file, sep='\t', index_col=0)
    
    # Convert counts to integers (required by PyDESeq2)
    counts_df = counts_df.round().astype(int)
    
    # remove low gene counts
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
    counts_df = counts_df[genes_to_keep]
    return counts_df, metadata

def get_comparisons(sample_info, custom_comparisons=None):
    """Get all pairwise comparisons or custom comparisons"""
    conditions = sample_info['Condition'].unique()
    print(f"Found conditions: {conditions}")
    
    if custom_comparisons:
        # Parse custom comparisons
        comparisons = []
        for comp in custom_comparisons.split(','):
            cond1, cond2 = comp.strip().split('_vs_')
            comparisons.append((cond1, cond2))
        print(f"Running custom comparisons: {comparisons}")
    else:
        # Generate all pairwise comparisons
        comparisons = list(combinations(conditions, 2))
        print(f"Running all pairwise comparisons: {comparisons}")
    
    return comparisons

def run_pydeseq2_base(counts_df, sample_info):
    """Run base PyDESeq2 analysis (preprocessing)"""
    print("Running PyDESeq2 base analysis...")
    
    # Create DeseqDataSet
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=sample_info,
        design_factors="Condition",
        refit_cooks=True
    )
    
    # Run the preprocessing
    print("Running DESeq2 preprocessing...")
    dds.deseq2()
    
    return dds

def run_comparison(dds, condition1, condition2):
    """Run a specific pairwise comparison"""
    print(f"Running comparison: {condition1} vs {condition2}")
    
    # Create contrast for this comparison
    contrast = ["Condition", condition1, condition2]
    
    # Get statistical results for this contrast
    stat_res = DeseqStats(dds, contrast=contrast)
    stat_res.summary()
    
    return stat_res

def create_results_dataframe(stat_res, comparison_name):
    """Create results dataframe for a specific comparison"""
    print(f"Creating results dataframe for {comparison_name}...")
    
    results_df = pd.DataFrame({
        'gene_id': stat_res.results_df.index,
        'baseMean': stat_res.results_df['baseMean'],
        'log2FoldChange': stat_res.results_df['log2FoldChange'],
        'lfcSE': stat_res.results_df['lfcSE'],
        'stat': stat_res.results_df['stat'],
        'pvalue': stat_res.results_df['pvalue'],
        'padj': stat_res.results_df['padj']
    })
    
    # Sort by adjusted p-value
    results_df = results_df.sort_values('padj')
    
    return results_df

def create_comparison_plots(dds, stat_res, results_df, comparison_name, plots_dir, padj_thresh, lfc_thresh):
    """Create plots for a specific comparison"""
    print(f"Creating plots for {comparison_name}...")
    
    # Create comparison-specific directory
    comp_plots_dir = f'{plots_dir}/{comparison_name}'
    comp_plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 1. MA Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Calculate mean expression and log2 fold change
    mean_expr = results_df['baseMean']
    log2fc = results_df['log2FoldChange']
    
    # Color by significance
    colors = ['red' if (padj < padj_thresh and abs(lfc) > lfc_thresh) 
              else 'gray' for padj, lfc in zip(results_df['padj'].fillna(1), log2fc)]
    
    ax.scatter(np.log10(mean_expr + 1), log2fc, c=colors, alpha=0.6, s=1)
    ax.axhline(y=lfc_thresh, color='blue', linestyle='--', alpha=0.7)
    ax.axhline(y=-lfc_thresh, color='blue', linestyle='--', alpha=0.7)
    ax.set_xlabel('Log10 Mean Expression')
    ax.set_ylabel('Log2 Fold Change')
    ax.set_title(f'MA Plot: {comparison_name}')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{comp_plots_dir}\ma_plot.pdf')
    plt.close()
    
    # 2. Volcano Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Prepare data
    log2fc = results_df['log2FoldChange']
    neg_log10_padj = -np.log10(results_df['padj'].fillna(1))
    
    # Color by significance
    significant = ((results_df['padj'] < padj_thresh) & 
                  (abs(results_df['log2FoldChange']) > lfc_thresh))
    
    colors = ['red' if sig else 'gray' for sig in significant]
    
    ax.scatter(log2fc, neg_log10_padj, c=colors, alpha=0.6, s=1)
    ax.axvline(x=lfc_thresh, color='blue', linestyle='--', alpha=0.7)
    ax.axvline(x=-lfc_thresh, color='blue', linestyle='--', alpha=0.7)
    ax.axhline(y=-np.log10(padj_thresh), color='blue', linestyle='--', alpha=0.7)
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10 Adjusted P-value')
    ax.set_title(f'Volcano Plot: {comparison_name}')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{comp_plots_dir}/volcano_plot.pdf')
    plt.close()

def create_global_plots(dds, sample_info, plots_dir):
    """Create plots that show all samples/conditions together"""
    print("Creating global plots...")
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Get normalized counts for plotting
    norm_counts = dds.layers["normed_counts"]
    
    # 1. PCA Plot
    print("Creating PCA plot...")
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Get log2 transformed data
    log2_counts = np.log2(norm_counts + 1)
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(log2_counts)
    
    # Create PCA dataframe
    pca_df = pd.DataFrame({
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1],
        'condition': sample_info['Condition'].values,
        'sample': sample_info.index
    })
    
    # Plot
    for condition in pca_df['condition'].unique():
        mask = pca_df['condition'] == condition
        ax.scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC2'], 
                  label=condition, s=100, alpha=0.7)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    ax.set_title('PCA of All Samples')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{plots_dir}/pca_all_samples.pdf')
    plt.close()
    
    # 2. Heatmap of top variable genes
    print("Creating heatmap...")
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Get top 50 most variable genes
    log2_counts_df = pd.DataFrame(log2_counts, index=dds.obs_names, columns=dds.var_names)
    gene_vars = log2_counts_df.var(axis=1)
    top_vars = gene_vars.nlargest(50).index
    
    # Create heatmap data
    heatmap_data = log2_counts_df.loc[top_vars]
    
    # Z-score normalization (row-wise)
    heatmap_data_zscore = heatmap_data.sub(heatmap_data.mean(axis=1), axis=0).div(heatmap_data.std(axis=1), axis=0)
    
    # Create color map for conditions
    unique_conditions = sample_info['Condition'].unique()
    condition_colors = {condition: color for condition, color in 
                       zip(unique_conditions, sns.color_palette("husl", len(unique_conditions)))}
    col_colors = [condition_colors[condition] for condition in sample_info['Condition']]
    
    # Plot heatmap
    sns.clustermap(heatmap_data_zscore, 
                   col_colors=col_colors,
                   cmap='RdBu_r', 
                   center=0,
                   figsize=(12, 8),
                   yticklabels=False)
    
    plt.suptitle('Heatmap of Top 50 Most Variable Genes (All Conditions)', y=1.02)
    plt.savefig(f'{plots_dir}\heatmap_top_genes_all.pdf', bbox_inches='tight')
    plt.close()

def create_comparison_summary(results_dict, outdir, padj_thresh, lfc_thresh):
    """Create summary for all comparisons"""
    print("Creating comparison summary...")
    
    summary_data = []
    for comparison_name, results_df in results_dict.items():
        total_genes = len(results_df)
        significant_genes = ((results_df['padj'] < padj_thresh) & 
                            (abs(results_df['log2FoldChange']) > lfc_thresh)).sum()
        upregulated = ((results_df['padj'] < padj_thresh) & 
                       (results_df['log2FoldChange'] > lfc_thresh)).sum()
        downregulated = ((results_df['padj'] < padj_thresh) & 
                         (results_df['log2FoldChange'] < -lfc_thresh)).sum()
        
        summary_data.append({
            'Comparison': comparison_name,
            'Total_Genes': total_genes,
            'Significant_Genes': significant_genes,
            'Upregulated': upregulated,
            'Downregulated': downregulated
        })
    
    # Create summary dataframe
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(f"{outdir}/comparison_summary.csv", index=False)
    
    # Create detailed text summary
    summary_text = f"""PyDESeq2 Multiple Conditions Analysis Summary
==========================================

Analysis Parameters:
- Adjusted p-value threshold: {padj_thresh}
- Log2 fold change threshold: {lfc_thresh}

Comparisons Performed:
"""
    
    for _, row in summary_df.iterrows():
        summary_text += f"""
{row['Comparison']}:
  - Total genes: {row['Total_Genes']}
  - Significant genes: {row['Significant_Genes']}
  - Upregulated: {row['Upregulated']}
  - Downregulated: {row['Downregulated']}
"""
    
    summary_text += f"""
Files Generated:
- comparison_summary.csv: Summary statistics for all comparisons
- For each comparison (e.g., condition1_vs_condition2/):
  - differential_expression_results.csv: Complete results
  - ma_plot.pdf: MA plot
  - volcano_plot.pdf: Volcano plot
- Global files:
  - normalized_counts.csv: Normalized counts for all samples
  - log2_transformed_counts.csv: Log2-transformed counts
  - pca_all_samples.pdf: PCA plot of all samples
  - heatmap_top_genes_all.pdf: Heatmap of top variable genes
"""
    
    with open(f'{outdir}/pydeseq2_summary.txt', 'w') as f:
        f.write(summary_text)

def main():
    args = parse_arguments()
    
    # Create output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    counts_df, sample_info = load_data(args.counts, args.samples)
    
    # Get comparisons to perform
    comparisons = get_comparisons(sample_info, args.comparisons)
    
    # Run base PyDESeq2 analysis (preprocessing)
    dds = run_pydeseq2_base(counts_df, sample_info)
    
    # Save global outputs
    print("Saving global outputs...")
    
    # Save normalized counts
    
    normalized_counts = pd.DataFrame(dds.layers["normed_counts"], 
                                       index=dds.obs_names, 
                                       columns=dds.var_names)
    
    normalized_counts.to_csv(f'{outdir}/normalized_counts.csv')
    
    # Save log2 transformed counts
    log2_counts = np.log2(normalized_counts + 1)
    log2_counts.to_csv(f'{outdir}/log2_transformed_counts.csv')
    
    # Create global plots
    create_global_plots(dds, sample_info, f'{outdir}/plots')
    
    # Run all comparisons
    results_dict = {}
    
    for condition1, condition2 in comparisons:
        comparison_name = f"{condition1}_vs_{condition2}"
        print(f"\n=== Processing {comparison_name} ===")
        
        # Run comparison
        stat_res = run_comparison(dds, condition1, condition2)
        
        # Create results dataframe
        results_df = create_results_dataframe(stat_res, comparison_name)
        
        # Store results
        results_dict[comparison_name] = results_df
        
        # Create comparison-specific directory
        comp_outdir = f'{outdir}/comparison_name'
        comp_outdir.mkdir(parents=True, exist_ok=True)
        
        # Save results for this comparison
        results_df.to_csv(f'{comp_outdir}/differential_expression_results.csv', index=False)
        
        # Create plots for this comparison
        create_comparison_plots(dds, stat_res, results_df, comparison_name, 
                              f'{outdir}/plots', args.padj, args.lfc)
    
    # Create overall summary
    create_comparison_summary(results_dict, outdir, args.padj, args.lfc)
    
    print("\nPyDESeq2 multiple conditions analysis completed successfully!")
    print(f"Results saved in: {outdir}")

if __name__ == "__main__":
    main()