#!/bin/bash

# hg19 
# This script pulls the latest RefSeq GFF annotation from NCBI and parses it into a format compatible with 
# the epitope ID scripts. Move the resulting files into the appropriate /pwd/annotation/ when complete

wget -N ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
gunzip -f GRCh37_latest_genomic.gff.gz
perl parsers/parse_NCBI_RefSeq_GFF.pl GRCh37_latest_genomic.gff temp.gff

echo "Tiling Genome..."
perl genome_bin/bin_genome.pl ~/Desktop/GENOMES/hg19/hg19.fa 500 hg19_BIN.gff
echo "Complete"

echo "Intersecting gene annotations..."
bedtools subtract -a hg19_BIN.gff -b temp.gff > hg19_BIN_temp.gff
echo "Complete"

echo "Merging annotation and bin file..."
perl genome_bin/rename_BIN_GFF.pl hg19_BIN_temp.gff hg19_BIN_filter.gff
cat hg19_BIN_filter.gff temp.gff > hg19_ALL.gff
sort -k 1,1 -k4,4n hg19_ALL.gff > genome_annotation.gff
gzip -f genome_annotation.gff
echo "Complete"

rm GRCh37_latest_genomic.gff temp.gff hg19_BIN.gff hg19_BIN_temp.gff hg19_BIN_filter.gff hg19_ALL.gff
