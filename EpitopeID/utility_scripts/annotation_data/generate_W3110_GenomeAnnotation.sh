#!/bin/bash

# W3110

# Required software:
# bedtools
# wget
# Perl 5.18+

usage()
{
    echo 'generate_W3110_GenomeAnnotation.sh -g /path/to/GenomeFASTA -b BIN size (bp)'
    echo 'eg: bash generate_W3110_GenomeAnnotation.sh -g /genome/W3110.fa -b 250'
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

perl parsers/parse_W3110_ORF_GFF.pl ORF/ec2_W3110_annotation.gff $BIN W3110_ORF.gff

echo "Tiling Genome..."
perl genome_bin/bin_genome.pl $GENOME $BIN W3110_BIN.gff
echo "Complete"

echo "Intersecting gene annotations..."
bedtools subtract -a W3110_BIN.gff -b W3110_ORF.gff > W3110_BIN_temp.gff
echo "Complete"

echo "Merging annotation and bin file..."
perl genome_bin/rename_BIN_GFF.pl W3110_BIN_temp.gff W3110_BIN_filter.gff
cat W3110_BIN_filter.gff W3110_ORF.gff > W3110_ALL.gff
sort -k 1,1 -k4,4n W3110_ALL.gff > genome_annotation.gff
gzip -f genome_annotation.gff
echo "Complete"

rm W3110_ORF.gff W3110_BIN.gff W3110_BIN_temp.gff W3110_BIN_filter.gff W3110_ALL.gff
