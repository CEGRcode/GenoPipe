#!/bin/bash

# This script simulates by sampling paired-end reads from a synthetic genome as BED coordinate
# using custom scrripts and then extracting the sequence using Bedtools. Then the FASTA has 
# quality scores added and the resulting FASTQ is compressed to fastq.gz files.

# Simulation index is like a replicate index and is used in the final output filename
# Depth specifies how many paired end reads to generate for the dataset
# Seed can be be set for the sampling step to generate the BED files
# Synthetic Genome FASTA file needs to be indicated (with inserted epitope, with deletion 
#		removed, with SNPs sprinkled into genome)
# Output FASTQ files are saved to /path/to/results/FASTQ/Simulation_X_R*.fastq.gz

# Required software:
# BWA v0.7.15
# bedtools v2.27.1
# samtools v1.5
# python 2.7.14
# pysam

usage()
{
	echo 'bash simulate.sh -i <simulationIndex> -d [10K|50K|100K|200K|500K|1M|2M|5M|7M|10M|50M] -s <seed> -g <syntheticGenomeFASTA> -o <outputDirectory>'
	echo 'eg: bash simulate.sh -i 1 -d 100K -s 1000 -g /path/to/syntheticgenome -o /path/to/results'
	exit
}

if [ "$#" -ne 10 ]; then
	usage
fi

declare -A GETDEPTH=( 
		["10K"]=10000 \
		["50K"]=50000 \
		["100K"]=100000 \
		["200K"]=200000 \
		["500K"]=500000 \
		["600K"]=600000 \
		["700K"]=700000 \
		["800K"]=800000 \
		["900K"]=900000 \
		["1M"]=1000000 \
		["1.5M"]=1500000 \
		["2M"]=2000000 \
		["3M"]=3000000 \
		["4M"]=4000000 \
		["5M"]=5000000 \
		["7M"]=7000000 \
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

RAND=../scripts/generate_random_BED_from_Genome_Size.pl
CONVERT=../scripts/convert_FASTA_to_GZIP-FASTQ.pl


perl $RAND $GENOME.fai $DEPTH $SEED $BED
bedtools getfasta -name -s -fi $GENOME -bed $BED\_R1.bed.gz -fo $FQ\_R1.fa
bedtools getfasta -name -s -fi $GENOME -bed $BED\_R2.bed.gz -fo $FQ\_R2.fa
perl $CONVERT $FQ\_R1.fa $FQ\_R1.fastq.gz
perl $CONVERT $FQ\_R2.fa $FQ\_R2.fastq.gz
rm ${FQ}_*.fa
