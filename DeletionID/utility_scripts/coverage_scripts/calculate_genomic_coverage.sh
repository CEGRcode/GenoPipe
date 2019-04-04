#!/bin/bash

# This script calculates the mappability of a user-specified series of genomic coordinates at a specified read
# length. The resulting mappability files should be moved into the appropriate /pwd/mappability/ folder
 
# Required software:
# BWA v0.7.14+
# bedtools v2.26+
# samtools v1.7+
# python 2.7.14
# pysam

usage()
{
    echo 'calculate_genomic_coverage.sh -f <Genomic FASTA file> -c <BED coordinate file> -r <Read Length> -t <Threads - Default 1>'
    echo 'eg: bash calculate_genomic_coverage.sh -f sacCer3.fa -c sacCer3_ORF.bed -r 50 -t 2'
    exit
}

if [ "$#" -ne 6 ] && [ "$#" -ne 8 ]; then
    usage
fi

THREAD=1

while getopts ":f:c:r:t:" IN; do
    case "${IN}" in
        f)
            GENOME=${OPTARG}
            ;;
        c)
            COORD=${OPTARG}
            ;;
        r)
            READLENGTH=${OPTARG}
            ;;
        t)
            THREAD=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${GENOME}" ] || [ -z "${COORD}" ] || [ -z "${READLENGTH}" ]; then
    usage
fi

echo "Genome = ${GENOME}"
echo "BED file = ${COORD}"
echo "Read length = ${READLENGTH}"
echo "Threads = ${THREAD}"

# Hard-coded relative paths
BORDER=coverage_scripts/expand_BED_by_Border.pl
TILE=coverage_scripts/tile_genome.pl 
COVERAGE=coverage_scripts/calculate_region_coverage.py

# Add pads to coordinates under consideration
perl $BORDER $COORD $READLENGTH temp.bed

# Get genomic FASTA sequence for padded coordinates
bedtools getfasta -fi $GENOME -bed temp.bed -fo temp.fa

# Tile genome at user-specified read length
perl $TILE temp.fa $READLENGTH TILE_GENOME.fa

# Align to reference genome, filtering for unique reads
bwa mem -t $THREAD $GENOME TILE_GENOME.fa | samtools view -Shb -q 5 - | samtools sort -o ALIGN_GENOME.bam -
samtools index ALIGN_GENOME.bam

# Calculate % uniquely mappable given user-specifed genomic coordinate file
python2 $COVERAGE -b ALIGN_GENOME.bam -c $COORD -o $READLENGTH\bp_Cov.out

# Clean up intermediary files
rm TILE_GENOME.fa ALIGN_GENOME.bam* temp.bed temp.fa 
