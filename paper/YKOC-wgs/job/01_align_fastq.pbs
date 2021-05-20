#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=16gb
#PBS -l walltime=01:00:00
#PBS -A open
#PBS -o logs/align-picard.log.out
#PBS -e logs/align-picard.log.err
#PBS -t 1-9010

module load gcc/8.3.1
module load samtools/1.5
module load bwa/0.7.15
module load picard/2.20.8 

WRK=/path/to/GenoPipe/paper/YKOC-wgs
cd $WRK

[ -d logs ] || mkdir logs
[ -d results/BAM ] || mkdir results/BAM

GENOME=../input/sacCer3.fa
METADATA=210403_PRJEB27160_accessions.txt
INDEX=$(($PBS_ARRAYID+1))

INFO=`sed "${INDEX}q;d" $METADATA`
ERR=`echo $INFO | awk '{print $4}'`
ERS=`echo $INFO | awk '{print $2}'`
FQ1=results/FASTQ/$ERR/$ERR\_1.fastq.gz
FQ2=results/FASTQ/$ERR/$ERR\_2.fastq.gz
BAM=$WRK/results/BAM/$ERS
#echo $INFO

start=`date +%s`
# Align with Pugh Lab standard alignment pipeline parameters
# -T INT  Don't output alignments with score lower than INT. This option affects output and occsaionally SAM flag 2.
# -h INT  If query has not more than INT hits with score higher than 880% of the best hit, output them all in the XA tag.
# -M      Mark shorter split hits as secondary (for Picard compatibility).
echo "(${PBS_ARRAYID}) Aligning $ERR fastq files..."
bwa mem -v 1 -T '30' -h '5' -t 4 -M $GENOME $FQ1 $FQ2 > $BAM.unsorted.bam
# Sort
echo "(${PBS_ARRAYID}) Sorting $ERR ..."
samtools sort $BAM.unsorted.bam > $BAM.unmarked.bam
# Mark duplicates
echo "(${PBS_ARRAYID}) Marking duplicates $ERR ..."
picard MarkDuplicates \
	INPUT=$BAM.unmarked.bam \
	OUTPUT=$BAM.marked.bam \
	METRICS_FILE=$BAM.metrics.txt \
	REMOVE_DUPLICATES='false' ASSUME_SORTED='true' \
	DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES' \
	#READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' \
	#OPTICAL_DUPLICATE_PIXEL_DISTANCE='100' \
	#VALIDATION_STRINGENCY='LENIENT' #VERBOSITY=ERROR
# Filter duplicates
echo "(${PBS_ARRAYID}) Filter duplicates $ERR ..."
samtools view  -h -b -f 0x1 -F 0x404 \
	-o $BAM.bam \
	$BAM.marked.bam
# Index
echo "(${PBS_ARRAYID}) Index $ERR ..."
samtools index $BAM.bam
end=`date +%s`
runtime=$((end-start))
echo "...alignment and indexing for ($PBS_ARRAYID) $ERR.fq > $ERS.bam finished in ${runtime}"
# Clean-up
rm $BAM.unsorted.bam $BAM.unmarked.bam $BAM.marked.bam
