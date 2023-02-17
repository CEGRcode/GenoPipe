#!/bin/bash

set -exo

# FIRST CHANGE PATH TO EXECUTE
WRK=/path/to/GenoPipe/paper/SyntheticEpitope
cd $WRK

# CALCULATE=scripts/calculate_detection_Stats.pl
TALLY=scripts/analyze_eid_results.py
BAR=scripts/build_barplots.py
LINE=scripts/build_lineplots.py

REF="sacCer3"
DEPTH="1M"
LENGTH="500"
EPITOPE=R$LENGTH
AFACTOR=Rap1
BFACTOR=Reb1
SUMMARY=results/MixSummaryReport_$REF\_$DEPTH.txt
RSUMMARY=results/MixRuntimeSummaryReport_$REF\_$DEPTH.txt
[ -f $SUMMARY ] && rm $SUMMARY
[ -f $RSUMMARY ] && rm $RSUMMARY

EPITOPE=R$LENGTH
for PER in 90 80 70 60 50 40 30 20 10
do
	CATRAW=results/mix_yeast/mix_yeast_$DEPTH\_$PER\_$AFACTOR\_$BFACTOR-all-IDs.tab
	head -n 9999 results/mix_yeast/$DEPTH/$PER\_$AFACTOR\_$BFACTOR/ID/*.tab > $CATRAW
	# perl $CALCULATE $CATRAW $PROTEIN RANDOM_SEQ_$LENGTH >> $SUMMARY/$DEPTH\_$EPITOPE\_$PROTEIN\_summary.txt
	python $TALLY -i results/mix_yeast/$DEPTH/$PER\_$AFACTOR\_$BFACTOR/ID/ -e RANDOM_SEQ_$LENGTH -t $AFACTOR \
		| sort > A-temp
	python $TALLY -i results/mix_yeast/$DEPTH/$PER\_$AFACTOR\_$BFACTOR/ID/ -e RANDOM_SEQ_$LENGTH -t $BFACTOR \
		| sort > B-temp
	paste A-temp B-temp >> $SUMMARY

	RUNTIMES=results/mix_yeast/mix_yeast_$DEPTH\_$PER\_$AFACTOR\_$BFACTOR-all-runtimes.tab
	grep 'finished in' results/mix_yeast/$DEPTH/$PER\_$AFACTOR\_$BFACTOR/runtime/*.runtime > $RUNTIMES
	paste <(cut -d":" -f1 $RUNTIMES) <(awk '{FS=" "}{OFS="\t"}{print $2,$6}' $RUNTIMES) >> $RSUMMARY
done

# Build yeast figures
python $LINE -i $SUMMARY -o results/ID-Mix-tally_$REF\_$DEPTH.png

REF="hg19"
DEPTH="20M"
LENGTH="500"
EPITOPE=R$LENGTH
AFACTOR=CTCF
BFACTOR=POLR2H
SUMMARY=results/MixSummaryReport_$REF\_$DEPTH.txt
RSUMMARY=results/MixRuntimeSummaryReport_$REF\_$DEPTH.txt
[ -f $SUMMARY ] && rm $SUMMARY
[ -f $RSUMMARY ] && rm $RSUMMARY

EPITOPE=R$LENGTH
for PER in 90 80 70 60 50 40 30 20 10
do
	CATRAW=results/mix_human/mix_human_$DEPTH\_$PER\_$AFACTOR\_$BFACTOR-all-IDs.tab
	head -n 9999 results/mix_human/$DEPTH/$PER\_$AFACTOR\_$BFACTOR/ID/*.tab > $CATRAW
	# perl $CALCULATE $CATRAW $PROTEIN RANDOM_SEQ_$LENGTH >> $SUMMARY/$DEPTH\_$EPITOPE\_$PROTEIN\_summary.txt
	python $TALLY -i results/mix_human/$DEPTH/$PER\_$AFACTOR\_$BFACTOR/ID/ -e RANDOM_SEQ_$LENGTH -t $AFACTOR \
		| sort > A-temp
	python $TALLY -i results/mix_human/$DEPTH/$PER\_$AFACTOR\_$BFACTOR/ID/ -e RANDOM_SEQ_$LENGTH -t $BFACTOR \
		| sort > B-temp
	paste A-temp B-temp >> $SUMMARY

	RUNTIMES=results/mix_human/mix_human_$DEPTH\_$PER\_$AFACTOR\_$BFACTOR-all-runtimes.tab
	grep 'finished in' results/mix_human/$DEPTH/$PER\_$AFACTOR\_$BFACTOR/runtime/*.runtime > $RUNTIMES
	paste <(cut -d":" -f1 $RUNTIMES) <(awk '{FS=" "}{OFS="\t"}{print $2,$6}' $RUNTIMES) >> $RSUMMARY
done

# Build human figures
python $LINE -i $SUMMARY -o results/ID-Mix-tally_$REF\_$DEPTH.png
