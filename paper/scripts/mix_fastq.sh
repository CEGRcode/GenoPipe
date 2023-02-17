#!/bin/bash

# This script takes the first user-specified "a" pepercentage of reads
# (assumes randomized order of FASTQ) of some FASTQ fileA and combines it
# with the first b percentage of reads in some FASTQ fileB to create a
# mixed population output FASTQ file.

# Required software:
# gzip
# seqtk

module load anaconda3
source activate genopipe

usage()
{
	echo 'bash mix_fastq.sh -p <percentOfA> -A <fastqBaseA> -B <fastqBaseB> -o <outBase>'
	echo 'eg: bash mix_fastq.sh -p 90 -a results/Reb1-Cterm_100K/Simulation_1 -b results/Rap1-Nterm_100K/Simulation_1 -o results/mix_yeast/90_Reb1-Rap1_10K/Simulation_1'
	exit
}

if [ "$#" -ne 8 ]; then
	usage
fi

while getopts ":p:a:b:o:" IN; do
	case "${IN}" in
	a)
		AFQ=${OPTARG}
		;;
	b)
		BFQ=${OPTARG}
		;;
	p)
		PER=${OPTARG}
		;;
	o)
		OFQ=${OPTARG}
		;;
	*)
		usage
		;;
	esac
done

shift $((OPTIND-1))

# get percentages as float
A_NUM=`echo $PER | awk '{print $1/100}'`
B_NUM=`echo $A_NUM | awk '{print 1 - $1}'`

echo "FASTQ_A = ${AFQ} (${A_NUM})"
echo "FASTQ_B = ${BFQ} (${B_NUM})"
echo "OUTPUT = ${OFQ}"

# Check that the files haven't already been generated
#if [[ ! -f $OFQ\_readnames.txt && \
#	-f $OFQ\_R1.fastq.gz && \
#	-f $OFQ\_R2.fastq.gz ]];
#then
#	echo "Both R1 and R2 files have already been generated for ${OFQ}"
#	exit
#fi

# Cleanup part-generated FQs if indicator file present from a previous run
[ -f $OFQ\_readnames.txt ] && rm $OFQ\_R1.fastq*
[ -f $OFQ\_readnames.txt ] && rm $OFQ\_R2.fastq*

# Check FASTQ files from depth simulations exist
[ -f $AFQ\_R1.fastq.gz ] || echo "A_FASTQ_R1 does not exist: ${AFQ}_R1.fastq.gz"
[ -f $AFQ\_R2.fastq.gz ] || echo "A_FASTQ_R2 does not exist: ${AFQ}_R2.fastq.gz"
[ -f $BFQ\_R1.fastq.gz ] || echo "B_FASTQ_R1 does not exist: ${BFQ}_R1.fastq.gz"
[ -f $BFQ\_R2.fastq.gz ] || echo "B_FASTQ_R2 does not exist: ${BFQ}_R2.fastq.gz"

# Get num of lines to retrieve from FileA (repeat for FileB) ( wc*frac but rest is to keep it a multiple of 4)
A_LINES=`gzip -dc $AFQ\_R1.fastq.gz | wc -l | awk -v A=$A_NUM '{print int($1*A/4)*4}'`
B_LINES=`gzip -dc $BFQ\_R1.fastq.gz | wc -l | awk -v B=$B_NUM '{print int($1*B/4)*4}'`
echo "A NLINES = $A_LINES"
echo "B NLINES = $B_LINES"

# Get first x number of read lines
gzip -dc $AFQ\_R1.fastq.gz | head -n $A_LINES - > $OFQ\_R1.fastq
gzip -dc $BFQ\_R1.fastq.gz | head -n $B_LINES - >> $OFQ\_R1.fastq

# Get FASTQ readnames from R1 file
grep '^@' $OFQ\_R1.fastq | sed 's/^@//' > $OFQ\_readnames.txt

# Get mate-pairs from respective R2 files
seqtk subseq $AFQ\_R2.fastq.gz $OFQ\_readnames.txt > $OFQ\_R2.fastq
seqtk subseq $BFQ\_R2.fastq.gz $OFQ\_readnames.txt >> $OFQ\_R2.fastq

# Zip and clean-up files
gzip -f $OFQ\_R1.fastq
gzip -f $OFQ\_R2.fastq

#rm $OFQ\_readnames.txt
