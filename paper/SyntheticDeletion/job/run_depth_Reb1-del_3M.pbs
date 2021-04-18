#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=16gb
#PBS -l walltime=00:30:00
#PBS -A open
#PBS -o logs/depth.Reb1.3M.log.out
#PBS -e logs/depth.Reb1.3M.log.err
#PBS -t 1-1000

# FIRST CHANGE PATH TO EXECUTE
WRK=/path/to/GenoPipe/paper/SyntheticDeletion
cd $WRK

module load gcc/8.3.1
module load bedtools/2.27.1
module load bwa/0.7.15
module load samtools/1.5
module load anaconda3
source activate genopipe

INFO=`sed "10q;d" depth_simulations.txt`
LOCUS=`awk '{print $1}'  <(echo $INFO)`
DEPTH=`awk '{print $2}'  <(echo $INFO)`
BASE=`awk '{print $3}'  <(echo $INFO)`

GENOME=synthetic_genome/$LOCUS\_del.fa
OUTPUT=results/$LOCUS\_$DEPTH
SEED=$(($BASE+$PBS_ARRAYID))

start=`date +%s`
bash ../scripts/simulate.sh -i $PBS_ARRAYID -d $DEPTH -s $SEED -g $GENOME -o $OUTPUT
bash ../scripts/align.sh -i $PBS_ARRAYID -g ../input/sacCer3.fa -o $OUTPUT -t 4
end=`date +%s`
runtime=$((end-start))
echo "${LOCUS} ${DEPTH} simulate in ${runtime}"
