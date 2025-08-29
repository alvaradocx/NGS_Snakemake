# NGS Snakemake Basic QC and Basic Analysis (In Progress)
This is a basic snakemake pipeline for processing fastq data from [Herpes Simplex Virus 1 Infection of Human Brain Organoids and Pancreatic Stem Cell-Islets Drives Organoid-Specific Transcripts Associated with Alzheimer’s Disease and Autoimmune Diseases](https://pmc.ncbi.nlm.nih.gov/articles/PMC11640215/). Fastq samples can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272361).
## Overview

This pipeline performs RNA-seq analysis including:
- Quality control (FastQC, MultiQC)
- Read trimming (fastp)
- Alignment to human + HSV-1 reference genomes (HISAT2)
- Read counting (HTSeq)
- Differential expression analysis (PyDESeq2) **IN PROGRESS**

## Requirements

- Conda/Mamba
- Minimum 16GB RAM
- ~450GB free disk space

## Installation

### Virtual Environment
```bash
conda env create -f environment.yaml
conda activate snakemake-pipeline
```

## Configuration

1. **Edit `config.yaml`** with your settings:
   - Input/output directories
   - Genome build (grch37 or grch38)
   - Email for NCBI downloads
   - Sample information
   - [Config Template](./snakemake/envs/config.yaml)

2. **Required files (User Provided):**
   - Place FASTQ files in `{input_dir}/data/{sample_id}/`
   - HSV-1 GTF annotation file
     - Can be downloaded from [NCBI Genome Tool](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027937515.1/)
   - SRA metadata CSV file (included in Fastq download from SRA)

## Usage

### Basic run:
```bash
snakemake --cores 8
```

### Dry run (recommended first):
```bash
snakemake --dry-run
```

### With conda environments:
```bash
snakemake --cores 8 --use-conda
```

## Pipeline Steps

1. **Quality Control**: FastQC on raw reads
2. **Trimming**: fastp for adapter removal and quality trimming
3. **Post-trim QC**: FastQC on trimmed reads
4. **Reference Preparation**: Download and combine human + HSV-1 genomes
5. **Alignment**: HISAT2 alignment to combined reference
6. **Post-alignment QC**: Qualimap analysis
7. **Read Counting**: HTSeq for gene expression quantification
8. **Differential Expression**: PyDESeq2 analysis with plots
9. **Summary Report**: MultiQC aggregated report

## Output Structure

```
out/
├── fastqc/          # Quality control reports
├── fastp/           # Trimmed reads
├── align/           # BAM alignment files
├── htseq/           # Gene count matrices
├── pydeseq2/        # DE analysis results and plots
└── multiqc/         # Summary report
```

## Key Output Files

- `multiqc/final_multiqc_report.html` - Overall QC summary
- `htseq/combined_counts_matrix.txt` - Gene expression matrix
- `pydeseq2/differential_expression_results.csv` - DE results
- `pydeseq2/plots/` - PCA, volcano, and other plots

## Sample Configuration

The pipeline expects samples organized as:
```
{input_dir}/data/
├── {sample1}/
│   ├── {id}_1.fastq.gz
│   └── {id}_2.fastq.gz
└── {sample2}/
    ├── {id}_1.fastq.gz
    └── {id}_2.fastq.gz
```

## Troubleshooting

**Common issues:**
- Ensure FASTQ files are properly named and in correct directories
- Check that genome build matches available references
- Verify sufficient disk space for large reference files
- Check conda environment activation

**Performance:**
- Adjust `--cores` based on available CPU
- Large reference downloads may take time on first run
- Consider using `--rerun-incomplete` for failed runs

## Citation

If you use this pipeline, please cite the original study:
[Herpes Simplex Virus 1 Infection of Human Brain Organoids...](https://pmc.ncbi.nlm.nih.gov/articles/PMC11640215/)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

Chelsea Alvarado
alvaradocx4@gmail.com
