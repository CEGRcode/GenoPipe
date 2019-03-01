#!/bin/bash

# This script calculates the mappability of a user-specified series of genomic coordinates at a specified read
# length. The resulting mappability files should be moved into the appropriate /pwd/mappability/ folder
 
# Required software:
# BWA v0.7.14+
# samtools v1.7+
# python 2.7.14
# pysam

usage()
{
    echo 'generate_sacCer3_Mappability_File.sh -f <Genomic FASTA file> -c <BED coordinate file> -t <Threads - Default 1>'
    echo 'eg: bash generate_sacCer3_Mappability_File.sh -f sacCer3.fa -c sacCer3_ORF.bed -t 2'
    exit
}

if [ "$#" -ne 6 ] && [ "$#" -ne 8 ]; then
    usage
fi

THREAD=1

while getopts ":f:c:t:" IN; do
    case "${IN}" in
        f)
            GENOME=${OPTARG}
            ;;
        c)
            COORD=${OPTARG}
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

if [ -z "${GENOME}" ] || [ -z "${COORD}" ]; then
    usage
fi

echo "Genome = ${GENOME}"
echo "BED file = ${COORD}"
echo "Threads = ${THREAD}"

# Hard-coded relative paths
COVERAGE=coverage_scripts/calculate_genomic_coverage.sh
MERGE=coverage_scripts/merge_Coverage_Files.sh

sh $COVERAGE -f $GENOME -c $COORD -r 36 -t $THREAD
sh $COVERAGE -f $GENOME -c $COORD -r 40 -t $THREAD
sh $COVERAGE -f $GENOME -c $COORD -r 50 -t $THREAD
sh $COVERAGE -f $GENOME -c $COORD -r 100 -t $THREAD

sh $MERGE

rm *Cov.out
