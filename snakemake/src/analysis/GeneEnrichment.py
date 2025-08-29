#!/usr/bin/env python3
"""
run_python_enrichment.py
Gene enrichment analysis using gseapy
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import gseapy as gp
import warnings
warnings.filterwarnings('ignore')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Gene enrichment analysis')
    parser.add_argument('--significant', required=True, help='All significant genes file')
    parser.add_argument('--upregulated', required=True, help='Upregulated genes file')
    parser.add_argument('--downregulated', required=True, help='Downregulated genes file')
    parser.add_argument('--outdir', required=True, help='Output directory')
    
    return parser.parse_args()

def convert_gene_ids(gene_list):
    """Convert gene IDs to gene symbols if needed"""
    # Remove version numbers from Ensembl IDs
    cleaned_genes = []
    for gene in gene_list:
        if isinstance(gene, str):
            # Remove version number if present
            cleaned_gene = gene.split('.')[0]
            cleaned_genes.append(cleaned_gene)
    
    return cleaned_genes

def run_go_enrichment(gene_list, output_dir, analysis_name):
    """Run GO enrichment analysis"""
    print(f"Running GO enrichment for {analysis_name}...")
    
    if len(gene_list) == 0:
        print(f"No genes provided for {analysis_name}")
        return None
    
    try:
        # Run GO enrichment
        go_results = gp.enrichr(gene_list=gene_list,
                               gene_sets=['GO_Biological_Process_2023',
                                         'GO_Molecular_Function_2023',
                                         'GO_Cellular_Component_2023'],
                               organism='Human',
                               outdir=f"{output_dir}/go_{analysis_name}",
                               cutoff=0.05)
        
        return go_results
    except Exception as e:
        print(f"GO enrichment failed for {analysis_name}: {e}")
        return None

def run_kegg_enrichment(gene_list, output_dir, analysis_name):
    """Run KEGG pathway enrichment"""
    print(f"Running KEGG enrichment for {analysis_name}...")
    
    if len(gene_list) == 0:
        print(f"No genes provided for {analysis_name}")
        return None
    
    try:
        # Run KEGG enrichment
        kegg_results = gp.enrichr(gene_list=gene_list,
                                 gene_sets=['KEGG_2021_Human'],
                                 organism='Human',
                                 outdir=f"{output_dir}/kegg_{analysis_name}",
                                 cutoff=0.05)
        
        return kegg_results
    except Exception as e:
        print(f"KEGG enrichment failed for {analysis_name}: {e}")
        return None

def combine_enrichment_results(enrichr_results, analysis_type, gene_set_type):
    """Combine enrichment results into a single dataframe"""
    combined_results = []
    
    if enrichr_results is None:
        return pd.DataFrame()
    
    for gene_set, results in enrichr_results.results.items():
        if len(results) > 0:
            results_df = results.copy()
            results_df['Gene_set'] = gene_set
            results_df['Analysis_type'] = analysis_type
            combined_results.append(results_df)
    
    if combined_results:
        return pd.concat(combined_results, ignore_index=True)
    else:
        return pd.DataFrame()

def create_enrichment_plots(go_results_all, kegg_results_all, plots_dir):
    """Create enrichment visualization plots"""
    print("Creating enrichment plots...")
    
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Gene Enrichment Analysis Results', fontsize=16)
    
    # Plot GO results
    if not go_results_all.empty:
        # GO Biological Process
        go_bp = go_results_all[go_results_all['Gene_set'].str.contains('Biological_Process')]
        if not go_bp.empty:
            top_go_bp = go_bp.nsmallest(10, 'Adjusted P-value')
            axes[0, 0].barh(range(len(top_go_bp)), -np.log10(top_go_bp['Adjusted P-value']))
            axes[0, 0].set_yticks(range(len(top_go_bp)))
            axes[0, 0].set_yticklabels([term[:50] + '...' if len(term) > 50 else term 
                                      for term in top_go_bp['Term']], fontsize=8)
            axes[0, 0].set_xlabel('-Log10(Adjusted P-value)')
            axes[0, 0].set_title('GO Biological Process')
            axes[0, 0].grid(True, alpha=0.3)
        
        # GO Molecular Function
        go_mf = go_results_all[go_results_all['Gene_set'].str.contains('Molecular_Function')]
        if not go_mf.empty:
            top_go_mf = go_mf.nsmallest(10, 'Adjusted P-value')
            axes[0, 1].barh(range(len(top_go_mf)), -np.log10(top_go_mf['Adjusted P-value']))
            axes[0, 1].set_yticks(range(len(top_go_mf)))
            axes[0, 1].set_yticklabels([term[:50] + '...' if len(term) > 50 else term 
                                      for term in top_go_mf['Term']], fontsize=8)
            axes[0, 1].set_xlabel('-Log10(Adjusted P-value)')
            axes[0, 1].set_title('GO Molecular Function')
            axes[0, 1].grid(True, alpha=0.3)
    
    # Plot KEGG results
    if not kegg_results_all.empty:
        top_kegg = kegg_results_all.nsmallest(15, 'Adjusted P-value')
        axes[1, 0].barh(range(len(top_kegg)), -np.log10(top_kegg['Adjusted P-value']))
        axes[1, 0].set_yticks(range(len(top_kegg)))
        axes[1, 0].set_yticklabels([term[:50] + '...' if len(term) > 50 else term 
                                  for term in top_kegg['Term']], fontsize=8)
        axes[1, 0].set_xlabel('-Log10(Adjusted P-value)')
        axes[1, 0].set_title('KEGG Pathways')
        axes[1, 0].grid(True, alpha=0.3)
    
    # Combined dot plot
    if not go_results_all.empty or not kegg_results_all.empty:
        combined_results = []
        if not go_results_all.empty:
            combined_results.append(go_results_all.nsmallest(5, 'Adjusted P-value'))
        if not kegg_results_all.empty:
            combined_results.append(kegg_results_all.nsmallest(5, 'Adjusted P-value'))
        
        if combined_results:
            combined_df = pd.concat(combined_results)
            top_combined = combined_df.nsmallest(15, 'Adjusted P-value')
            
            # Create dot plot
            scatter = axes[1, 1].scatter(top_combined['Odds Ratio'], 
                                       range(len(top_combined)),
                                       s=-np.log10(top_combined['Adjusted P-value']) * 20,
                                       c=-np.log10(top_combined['Adjusted P-value']),
                                       cmap='viridis', alpha=0.7)
            
            axes[1, 1].set_yticks(range(len(top_combined)))
            axes[1, 1].set_yticklabels([term[:40] + '...' if len(term) > 40 else term 
                                      for term in top_combined['Term']], fontsize=8)
            axes[1, 1].set_xlabel('Odds Ratio')
            axes[1, 1].set_title('Top Enriched Terms (Dot Plot)')
            axes[1, 1].grid(True, alpha=0.3)
            
            # Add colorbar
            plt.colorbar(scatter, ax=axes[1, 1], label='-Log10(Adjusted P-value)')
    
    plt.tight_layout()
    plt.savefig(f'{plots_dir}/enrichment_plots.pdf', bbox_inches='tight', dpi=300)
    plt.close()

def main():
    args = parse_arguments()
    
    # Create output directories
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plots_dir = f'{outdir}/plots'
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Load gene lists
    print("Loading gene lists...")
    significant_genes = pd.read_csv(args.significant)
    upregulated_genes = pd.read_csv(args.upregulated)
    downregulated_genes = pd.read_csv(args.downregulated)