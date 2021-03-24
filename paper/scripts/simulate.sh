#!/bin/bash

# This script simulates data and runs EpitopeID on the results. You can point to the epitope coordinates,
# sequencing depth, starting seed, and reference genome to use.

# Required software:
# BWA v0.7.15
# bedtools v2.27.1
# samtools v1.5
# python 2.7.14
# pysam

usage()
{
	echo 'bash simulate.sh -i <simulationIndex> -d [10K|100K|1M|10M|50M] -s <seed> -g <syntheticGenomeFASTA> -o <outputDirectory> -t <numThreads>'
	echo 'eg: bash simulate.sh -i 1 -d 100K -s 1000 -g /path/to/genomewepitope -o /path/to/results -t 4'
	exit
}

if [ "$#" -ne 10 ] && [ "$#" -ne 12 ]; then
	usage
fi

THREAD=4
declare -A GETDEPTH=( ["10K"]=10000 \
		["100K"]=100000 \
		["1M"]=1000000 \
		["10M"]=10000000 \
		["20M"]=20000000 \
		["50M"]=50000000  )

while getopts ":i:d:s:g:o:t:" IN; do
	case "${IN}" in
	i)
		INDEX=${OPTARG}
		;;
	d)
		S_DEPTH=${OPTARG}
		;;
	s)
		SEED=${OPTARG}
		;;
	g)
		GENOME=${OPTARG}
		;;
	o)
		OUTDIR=${OPTARG}
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

if [[ ! " ${!GETDEPTH[@]} " =~ " ${S_DEPTH} " ]]; then
	usage
fi
DEPTH=${GETDEPTH[$S_DEPTH]}

echo "Simulation-i = ${INDEX}"
echo "Depth = ${DEPTH}"
echo "SeedBaseVal = ${SEED}"
echo "SynthGenome = ${GENOME}"
echo "ResultsDirectory = ${OUTDIR}"

BED=$OUTDIR/BED/Simulation_$INDEX
FQ=$OUTDIR/FASTQ/Simulation_$INDEX

[ -d $OUTDIR ] || mkdir $OUTDIR
[ -d $OUTDIR/BED ]   || mkdir $OUTDIR/BED
[ -d $OUTDIR/FASTQ ] || mkdir $OUTDIR/FASTQ

RAND=../scripts/generate_random_BED_from_Genomic_FASTA.pl
CONVERT=../scripts/convert_FASTA_to_GZIP-FASTQ.pl

perl $RAND $GENOME $DEPTH $SEED $BED
bedtools getfasta -name -s -fi $GENOME -bed $BED\_R1.bed.gz -fo $FQ\_R1.fa
bedtools getfasta -name -s -fi $GENOME -bed $BED\_R2.bed.gz -fo $FQ\_R2.fa
perl $CONVERT $FQ\_R1.fa $FQ\_R1.fastq.gz
perl $CONVERT $FQ\_R2.fa $FQ\_R2.fastq.gz
rm ${FQ}*.fa
