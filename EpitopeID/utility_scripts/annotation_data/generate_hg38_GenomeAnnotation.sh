#!/bin/bash

# hg38
# This script pulls the latest RefSeq GFF annotation from NCBI and parses it into a format compatible with 
# the epitope ID scripts. Move the resulting files into the appropriate /pwd/annotation/ when complete

usage()
{
    echo 'generate_hg38_GenomeAnnotation.sh -g /path/to/GenomeFASTA -b BIN size (bp)'
    echo 'eg: bash generate_hg38_GenomeAnnotation.sh -g /genome/sacCer3.fa -b 250'
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

wget -N ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
gunzip -f GRCh38_latest_genomic.gff.gz
perl parsers/parse_NCBI_RefSeq_GFF.pl GRCh38_latest_genomic.gff $BIN temp.gff

echo "Tiling Genome..."
perl genome_bin/bin_genome.pl $GENOME $BIN hg38_BIN.gff
echo "Complete"

echo "Intersecting gene annotations..."
bedtools subtract -a hg38_BIN.gff -b temp.gff > hg38_BIN_temp.gff
echo "Complete"

echo "Merging annotation and bin file..."
perl genome_bin/rename_BIN_GFF.pl hg38_BIN_temp.gff hg38_BIN_filter.gff
cat hg38_BIN_filter.gff temp.gff > hg38_ALL.gff
sort -k 1,1 -k4,4n hg38_ALL.gff > genome_annotation.gff
gzip -f genome_annotation.gff
echo "Complete"

rm GRCh38_latest_genomic.gff temp.gff hg38_BIN.gff hg38_BIN_temp.gff hg38_BIN_filter.gff hg38_ALL.gff
