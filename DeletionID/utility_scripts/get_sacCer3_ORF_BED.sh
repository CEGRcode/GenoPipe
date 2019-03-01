#!/bin/bash
  
# This script retrieves the latest sacCer3 feature file from SGD and parses it into BED format for downstream
# processing. The resulting BED files should be moved into the appropriate /pwd/genomic_coord/ folder

# Required software:
# wget

# Retrieve latest SGD features file in GFF format
wget -c https://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff

# Parse file to BED format while filtering to only be Verified ORFs
PARSE=parsers/convert_SGD-GFF_to_BED.pl
perl $PARSE saccharomyces_cerevisiae.gff sacCer3_VerifiedORF.bed

# Clean up temporary files
rm saccharomyces_cerevisiae.gff
