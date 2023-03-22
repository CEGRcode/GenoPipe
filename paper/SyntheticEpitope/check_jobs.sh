#!/bin/bash

# FIRST CHANGE PATH TO EXECUTE
WRK=/path/to/GenoPipe/paper/SyntheticEpitope
cd $WRK

N_EXPERIMENTS=`wc -l depth_simulations.txt |awk '{print $1}'`

for EPITOPE in "R500" "R100" "R50" "R20";
do
	REF="hg19"
	for PROTEIN in "MED12" "YY1" "USF1";
	do
		for DEPTH in "50M";
		do
			for i in {1..100};
			do
				FQ=results/$REF/$DEPTH/$EPITOPE/$PROTEIN/FASTQ/Simulation_$i
				#FQ=results/hg19/50M/R500/MED12/FASTQ/Simulation_$i
				[ -f $FQ\_R1.fastq.gz ] || echo $FQ\_R1.fastq.gz
				[ -f $FQ\_R2.fastq.gz ] || echo $FQ\_R2.fastq.gz
			done
		done
	done

	REF="sacCer3"
	for PROTEIN in "Reb1" "Sua7" "Taf2" "Spt4" "Spt7" "Gcn5" "Hsf1";
	do
		for DEPTH in "10M";
		do
			for i in {1..1000};
			do
				FQ=results/$REF/$DEPTH/$EPITOPE/$PROTEIN/FASTQ/Simulation_$i
				#FQ=results/sacCer3/10M/R500/Reb1/FASTQ/Simulation_$i
				[ -f $FQ\_R1.fastq.gz ] || echo $FQ\_R1.fastq.gz
				[ -f $FQ\_R2.fastq.gz ] || echo $FQ\_R2.fastq.gz
			done
		done
	done
done
