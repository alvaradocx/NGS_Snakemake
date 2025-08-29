#!/usr/bin/env python3
import sys
import os
from Bio import Entrez, SeqIO

email = sys.argv[1]
outdir= sys.argv[2]

# email needed for entrez access
Entrez.email=email

# extract genbank record
handle=Entrez.efetch(db="nuccore", id="JQ673480.1",rettype="gb",retmode='text')
rec_text = SeqIO.read(handle,"genbank")

# write out file
with open(os.path.join(outdir,"Herpesvirus_1_strain_KOS.fasta"),"w") as f:
    SeqIO.write(rec_text,f,"fasta")
