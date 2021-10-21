#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=16gb
#PBS -l walltime=02:00:00
#PBS -A open
#PBS -o logs/depth.tally.log.out
#PBS -e logs/depth.tally.log.err

module load anaconda3
source activate genopipe

WRK=/path/to/GenoPipe/paper/SyntheticStrain
cd $WRK

TALLY=scripts/parse_simulation_results.py
RUNTIME=scripts/parse_runtimes.py
# Parse hg19 results
for STRAIN in "HELA" "K562";
do
	for DEPTH in "1M" "2M" "5M" "10M" "20M";
	do
		DIR=results/hg19_$STRAIN\_$DEPTH
		echo "Tally for $DIR..."
		python $TALLY -v ../db/hg19_VCF/ -i $DIR/ID > $DIR\_scores.txt
		python $RUNTIME -i <(grep 'real' $DIR/ID/*.time) > $DIR\_runtimes.txt
	done
done
# Parse sacCer3 results
for STRAIN in "CEN.PK2-1Ca" "RM11-1A";
do
	for DEPTH in "10K" "50K" "100K" "500K" "1M" "2M";
	do
		DIR=results/sacCer3_$STRAIN\_$DEPTH
		echo "Tally for $DIR..."
		python $TALLY -v ../db/sacCer3_VCF/ -i $DIR/ID > $DIR\_scores.txt
		python $RUNTIME -i <(grep 'real' $DIR/ID/*.time) > $DIR\_runtimes.txt
	done
done
