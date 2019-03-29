#!/bin/bash
  
# sacCer3 
# This script pulls the latest transcript annotations from SGD and parses it into a format compatible with 
# the epitope ID scripts. Move the resulting files into the appropriate /pwd/annotation/ when complete

# Required software:
# bedtools
# wget
# Perl 5.18+

usage()
{
    echo 'generate_sacCer3_GenomeAnnotation.sh -g /path/to/GenomeFASTA -b BIN size (bp)'
    echo 'eg: bash generate_sacCer3_GenomeAnnotation.sh -g /genome/sacCer3.fa -b 250'
    exit
}

if [ "$#" -ne 4 ]; then
    usage
fi

while getopts ":g:b:" IN; do
    case "${IN}" in
        g)
            GENOME=${OPTARG}
            ;;
	b)
	    BIN=${OPTARG}
	    ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${GENOME}" ] || [ -z "${BIN}" ]; then
    usage
fi

echo "Genome FASTA file = ${GENOME}"
echo "Genome BIN size = ${BIN}"

wget -N https://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff
perl parsers/parse_sacCer3_ORF_GFF.pl saccharomyces_cerevisiae.gff $BIN sacCer3_ORF.gff

echo "Tiling Genome..."
perl genome_bin/bin_genome.pl $GENOME $BIN sacCer3_BIN.gff
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
