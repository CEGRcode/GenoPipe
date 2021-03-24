#!/bin/bash

# This script aligns simulated data then sorts and indexs the resulting BAM file in the 
# appropriate simulation results directory.

# Required software:
# BWA v0.7.15
# bedtools v2.27.1
# samtools v1.5
# python 3.8

usage()
{
	echo 'bash align.sh -i <simulationIndex> -g <syntheticGenomeFASTA> -o <outputDirectory> -t <numThreads>'
	echo 'eg: bash align.sh -i 1 -g /path/to/genome -o /path/to/results -t 4'
	exit
}

if [ "$#" -ne 6 ] && [ "$#" -ne 8 ]; then
	usage
fi

THREAD=4

while getopts ":i:g:o:t:" IN; do
	case "${IN}" in
	i)
		INDEX=${OPTARG}
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

echo "Simulation_i = ${INDEX}"
echo "ReferenceGenome = ${GENOME}"
echo "ResultsDirectory = ${OUTDIR}"

FQ=$OUTDIR/FASTQ/Simulation_$INDEX
BAM=$OUTDIR/BAM/Simulation_$INDEX

[ -d $OUTDIR/BAM ]   || mkdir $OUTDIR/BAM

bwa mem -t $THREAD $GENOME $FQ\_R1.fastq.gz $FQ\_R2.fastq.gz \
	| samtools sort \
	> $BAM.bam
samtools index $BAM.bam
