#!/bin/bash
  
# This script concatenates a series of FASTA files containing the nucleotide sequences of various epitopes
# and creates a bwa index. The resulting files should be moved into the appropriate /pwd/tag_FASTA/ folder

# Required software:
# Bowtie2 v2.2.5+ 

# Remove existing master tag file and bwa indexes if they exist
rm -f ALL_TAG.fa*
# Concatenate various posible FASTA files into master index
cat *.fa *.fna *.ffn *.fasta > ALL_TAG.fa
# Bowtie2 index command
bowtie2-build ALL_TAG.fa ALL_TAG.fa
