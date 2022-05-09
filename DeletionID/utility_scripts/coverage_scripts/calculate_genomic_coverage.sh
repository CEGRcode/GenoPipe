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
    echo 'calculate_genomic_coverage.sh -f <Genomic FASTA file> -c <BED coordinate file> -r <Read Length> -s <Step resolution for tiling> -t <Threads - Default 1>'
    echo 'eg: bash calculate_genomic_coverage.sh -f sacCer3.fa -c sacCer3_ORF.bed -r 50 -s 25 -t 2'
    exit
}

if [ "$#" -ne 8 ] && [ "$#" -ne 10 ]; then
    usage
fi

THREAD=1

while getopts ":f:c:r:s:t:" IN; do
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
	s)
	    STEP=${OPTARG}
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

if [ -z "${GENOME}" ] || [ -z "${COORD}" ] || [ -z "${READLENGTH}" ] || [ -z "${STEP}" ]; then
    usage
fi

echo "Genome = ${GENOME}"
echo "BED file = ${COORD}"
echo "Read length = ${READLENGTH}"
echo "Tiling step (bp) = ${STEP}"
echo "Threads = ${THREAD}"

# Hard-coded relative paths
BORDER=coverage_scripts/expand_BED_by_Border.pl
TILE=coverage_scripts/tile_genome.pl 
COVERAGE=coverage_scripts/calculate_region_coverage.py

# Add pads to coordinates under consideration
perl $BORDER $COORD $READLENGTH coord-pad.bed

# Get genomic FASTA sequence for padded coordinates
bedtools getfasta -fi $GENOME -bed coord-pad.bed -fo coord-pad.fa

# Tile genome at user-specified read length
perl $TILE coord-pad.fa $READLENGTH $STEP TILE_GENOME.fa

# Align to reference genome, filtering for unique reads
bwa mem -t $THREAD $GENOME TILE_GENOME.fa | samtools view -Shb -q 5 - | samtools sort -o ALIGN_GENOME.bam -
samtools index ALIGN_GENOME.bam

# Calculate % uniquely mappable given user-specifed genomic coordinate file
python $COVERAGE -b ALIGN_GENOME.bam -c $COORD -r $STEP -o $READLENGTH\bp_Cov.out

# Clean up intermediary files
rm TILE_GENOME.fa ALIGN_GENOME.bam* coord-pad.*
