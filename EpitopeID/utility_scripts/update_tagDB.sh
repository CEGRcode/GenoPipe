#!/bin/bash
  
# This script concatenates a series of FASTA files containing the nucleotide sequences of various epitopes
# and creates a bwa index. The resulting files should be moved into the appropriate /pwd/tag_FASTA/ folder

# Required software:
# BWA v0.7.14+

# Remove existing master tag file and bwa indexes if they exist
rm -f ALL_TAG.fa*
# Concatenate various posible FASTA files into master index
cat *.fa *.fna *.ffn *.fasta > ALL_TAG.fa
# BWA index command
bwa index ALL_TAG.fa


