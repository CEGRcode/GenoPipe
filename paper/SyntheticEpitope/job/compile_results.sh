#!/bin/bash

set -exo

# FIRST CHANGE PATH TO EXECUTE
WRK=/path/to/GenoPipe/paper/SyntheticEpitope
cd $WRK

# CALCULATE=scripts/calculate_detection_Stats.pl
TALLY=scripts/analyze_eid_results.py
BAR=scripts/build_barplots.py
VIOLIN=scripts/build_violinplots.py


REF="hg19"
SUMMARY=results/SummaryReport_$REF.txt
RSUMMARY=results/RuntimeSummaryReport_$REF.txt
[ -f $SUMMARY ] && rm $SUMMARY
[ -f $RSUMMARY ] && rm $RSUMMARY

for LENGTH in "500" "100" "50" "20";
do
	EPITOPE=R$LENGTH
	for PROTEIN in "CTCF" "POLR2H" "YY1";
	do
		for DEPTH in "100K" "1M" "10M" "20M" "50M";
		do
			CATRAW=results/$REF/$REF\_$DEPTH\_$EPITOPE\_$PROTEIN-all-IDs.tab
			head -n 9999 results/$REF/$DEPTH/$EPITOPE/$PROTEIN/ID/*.tab > $CATRAW
			# perl $CALCULATE $CATRAW $PROTEIN RANDOM_SEQ_$LENGTH >> $SUMMARY/$DEPTH\_$EPITOPE\_$PROTEIN\_summary.txt
			python $TALLY -i results/$REF/$DEPTH/$EPITOPE/$PROTEIN/ID/ -e RANDOM_SEQ_$LENGTH -t $PROTEIN >> $SUMMARY

			RUNTIMES=results/$REF/$REF\_$DEPTH\_$EPITOPE\_$PROTEIN-all-runtimes.tab
			grep 'finished in' results/$REF/$DEPTH/$EPITOPE/$PROTEIN/runtime/*.runtime > $RUNTIMES
			paste <(cut -d":" -f1 $RUNTIMES) <(awk '{FS=" "}{OFS="\t"}{print $2,$6}' $RUNTIMES) >> $RSUMMARY
		done
	done
done

# Build human figures
python $BAR -i $SUMMARY -o results/ID-tally_$REF.png
python $VIOLIN -i $RSUMMARY -o results/Runtimes_$REF.png

REF="sacCer3"
SUMMARY=results/SummaryReport_$REF.txt
RSUMMARY=results/RuntimeSummaryReport_$REF.txt
[ -f $SUMMARY ] && rm $SUMMARY
[ -f $RSUMMARY ] && rm $RSUMMARY

for LENGTH in "500" "100" "50" "20";
do
	EPITOPE=R$LENGTH
	#for PROTEIN in "Reb1" "Rap1" "Sua7" "Taf2" "Spt4" "Spt7" "Gcn5" "Hsf1" "Fzo1" "Lge1";
	for PROTEIN in "Reb1" "Rap1" "Sua7";
	do
		for DEPTH in "10K" "100K" "1M" "10M";
		do
			CATRAW=results/$REF/$REF\_$DEPTH\_$EPITOPE\_$PROTEIN-all-IDs.tab
			head -n 9999 results/$REF/$DEPTH/$EPITOPE/$PROTEIN/ID/*.tab > $CATRAW
			# perl $CALCULATE $CATRAW $PROTEIN RANDOM_SEQ_$LENGTH >> $SUMMARY/$DEPTH\_$EPITOPE\_$PROTEIN\_summary.txt
			python $TALLY -i results/$REF/$DEPTH/$EPITOPE/$PROTEIN/ID/ -e RANDOM_SEQ_$LENGTH -t $PROTEIN >> $SUMMARY

			RUNTIMES=results/$REF/$REF\_$DEPTH\_$EPITOPE\_$PROTEIN-all-runtimes.tab
			grep 'finished in' results/$REF/$DEPTH/$EPITOPE/$PROTEIN/runtime/*.runtime > $RUNTIMES
			paste <(cut -d":" -f1 $RUNTIMES) <(awk '{FS=" "}{OFS="\t"}{print $2,$6}' $RUNTIMES) >> $RSUMMARY
		done
	done
done

# Build yeast figures
python $BAR -i $SUMMARY -y -o results/ID-tally_$REF.png
python $VIOLIN -i $RSUMMARY -y -o results/Runtimes_$REF.png
