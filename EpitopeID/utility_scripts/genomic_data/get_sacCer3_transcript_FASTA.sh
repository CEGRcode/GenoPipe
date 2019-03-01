#!/bin/bash

# sacCer3 
# This script pulls the latest transcript FASTA from SGD and parses it into a format compatible with 
# the epitope ID scripts. Move the resulting files into the appropriate /pwd/gene_FASTA/ when complete

# Required software:
# wget
# BWA v0.7.14+

wget https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
gunzip -f orf_coding.fasta.gz
perl parsers/parse_sacCer3_ORF_FASTA.pl orf_coding.fasta gene.fasta
bwa index gene.fasta
rm orf_coding.fasta
