#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=16gb
#PBS -l walltime=02:00:00
#PBS -A open
#PBS -o logs/align.data.log.out
#PBS -e logs/align.data.log.err
#PBS -t 1-10

module load gcc
module load samtools
module load bwa
module load anaconda3
source activate ~/work/myconda/genopipe/

# FIRST CHANGE PATH TO EXECUTE
WRK=/path/to/GenoPipe/paper/BY4742-chipseq
cd $WRK

[ -d logs ] || mkdir logs
[ -d results/BAM ] || mkdir -p results/BAM
[ -d results/uniq-BAM ] || mkdir -p results/uniq-BAM

YGENOME=$WRK/../input/sacCer3.fa
CSGENOME=$WRK/../input/sacCer3_index

INDEX=$(($PBS_ARRAYID+1))

METADATA=SraRunInfo.csv
INFO=`sed "${INDEX}q;d" $METADATA`
SRR=`echo $INFO | cut -d"," -f1`
SAMPLE=`echo $INFO | cut -d"," -f12`
PLATFORM=`echo $INFO | cut -d"," -f19`
#PAIR=`echo $INFO | cut -d"," -f16`
#echo $INFO

FQ=$WRK/results/FASTQ/$SRR
BAM=$WRK/results/BAM/$SAMPLE

echo "($PBS_ARRAYID) Aligned $SRR $PLATFORM reads > $BAM"
if [[ " $PLATFORM " =~ " ABI_SOLID " ]]; then
	bowtie -C -S $CSGENOME <(gzip -dc $YGENOME $FQ\_1.fastq.gz) \
		| samtools sort \
		> $BAM.bam
	echo "(PBS_ARRAYID) $BAM single aligned (bowtie color space)"
elif [[ " $PLATFORM " =~ " ILLUMINA " ]]; then
	bwa mem $YGENOME $FQ\_1.fastq.gz $FQ\_2.fastq.gz -t 4 \
		| samtools sort \
		> $BAM.bam
	echo "($PBS_ARRAYID) $BAM pair aligned (BWA)"
fi

#samtools view -b -F4 $BAM > $WRK/results/uniq-BAM/$SAMPLE.bam

echo "($PBS_ARRAYID) Indexing..."
samtools index $BAM.bam
echo "($PBS_ARRAYID) Complete!"
