#!/usr/bin/bash

cd ..

#!/usr/bin/bash

ENSEMBL_RELEASE=75
ENSEMBL_GRCh37_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna
F=Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
ENSEMBL_GRCh37=${ENSEMBL_GRCh37_BASE}/${F}
OUT="/home/jupyter/snakemake/src/data/reference"

# get grCh37 from ensembl
wget $ENSEMBL_GRCh37.gz

# unzip file
gzip -d $F.gz

mv $F $OUT/grch37_genome.fa

rm $F