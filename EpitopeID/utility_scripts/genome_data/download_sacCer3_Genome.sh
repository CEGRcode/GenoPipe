#!/bin/bash
  
# sacCer3 
# This script downloads the sacCer3 genome from SGD and indexes it for use with BWA.
# Move the resulting files into the appropriate /pwd/genome/ when complete

# Required software:
# wget
# Perl 5.18+
# bwa v0.7.14+

wget -N https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz
tar -xvzf S288C_reference_genome_R64-1-1_20110203.tgz

echo "Parsing genome..."
perl parsers/parse_sacCer3_Genome_FASTA.pl S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa genome.fa
echo "Complete"

echo "BWA Indexing genome..."
bwa index genome.fa
echo "Complete"

rm S288C_reference_genome_R64-1-1_20110203.tgz
rm -r S288C_reference_genome_R64-1-1_20110203/
