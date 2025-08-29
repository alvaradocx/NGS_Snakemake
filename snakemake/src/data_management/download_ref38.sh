#!/usr/bin/bash

ENSEMBL_RELEASE=84
ENSEMBL_GRCh38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna
F=Homo_sapiens.GRCh38.dna.primary_assembly.fa
ENSEMBL_GRCh38=${ENSEMBL_GRCh38_BASE}/${F}
OUT="/home/jupyter/snakemake/src/data/reference"

# get grCh38 from ensembl
wget $ENSEMBL_GRCh38.gz

# unzip file
gzip -d $F.gz

mv $F $OUT/grch38_genome.fa

rm $F
