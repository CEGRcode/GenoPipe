#!/bin/bash
  
# sacCer3 
# This script pulls the latest transcript annotations from SGD and parses it into a format compatible with 
# the epitope ID scripts. Move the resulting files into the appropriate /pwd/annotation/ when complete

# Required software:
# bedtools
# wget
# Perl 5.18+

wget -N https://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff
perl parsers/parse_sacCer3_ORF_GFF.pl saccharomyces_cerevisiae.gff sacCer3_ORF.gff

echo "Tiling Genome..."
perl genome_bin/bin_genome.pl ~/Desktop/GENOMES/sacCer3/sacCer3.fa 500 sacCer3_BIN.gff
echo "Complete"

echo "Intersecting gene annotations..."
bedtools subtract -a sacCer3_BIN.gff -b sacCer3_ORF.gff > sacCer3_BIN_temp.gff
echo "Complete"

echo "Merging annotation and bin file..."
perl genome_bin/rename_BIN_GFF.pl sacCer3_BIN_temp.gff sacCer3_BIN_filter.gff
cat sacCer3_BIN_filter.gff sacCer3_ORF.gff > sacCer3_ALL.gff
sort -k 1,1 -k4,4n sacCer3_ALL.gff > genome_annotation.gff
gzip -f genome_annotation.gff
echo "Complete"

rm saccharomyces_cerevisiae.gff sacCer3_ORF.gff sacCer3_BIN.gff sacCer3_BIN_temp.gff sacCer3_BIN_filter.gff sacCer3_ALL.gff
