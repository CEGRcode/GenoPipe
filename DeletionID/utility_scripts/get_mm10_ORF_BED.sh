#!/bin/bash
  
# This script retrieves the latest mm10 feature file from UCSC and parses it into BED format for downstream
# processing. The resulting BED files should be moved into the appropriate /pwd/mm10_Del/ folder

# Required software:
# wget
# Perl

# Retrieve latest SGD features file in GFF format
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gff.gz

# Parse file to BED format while filtering
perl parsers/parse_NCBI_RefSeq_GFF.pl GCF_000001635.26_GRCm38.p6_genomic.gff.gz coord.bed

# Clean up temporary files
rm GCF_000001635.26_GRCm38.p6_genomic.gff.gz
