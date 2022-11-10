#!/bin/bash

set -exo

# FIRST CHANGE PATH TO EXECUTE
WRK=/path/to/GenoPipe/paper/SyntheticEpitope
cd $WRK

# CALCULATE=scripts/calculate_detection_Stats.pl
TALLY=scripts/analyze_eid_results.py

# SUMMARY=results/sacCer3/SummaryReports
# [ -d $SUMMARY ] || mkdir $SUMMARY
SUMMARY=results/sacCer3/SummaryReport.txt
[ -f $SUMMARY ] && rm $SUMMARY

for LENGTH in "500" "100" "50" "20";
do
        EPITOPE=R$LENGTH
        for PROTEIN in "Reb1" "Rap1" "Sua7" "Taf2" "Spt4" "Spt7" "Gcn5" "Hsf1" "Fzo1" "Lge1";
        do
                #break
                REF="sacCer3"
                for DEPTH in "10K" "100K" "1M" "10M";
                do
                      CATRAW=results/$REF/$DEPTH/$EPITOPE/$PROTEIN/$REF\_$DEPTH\_$EPITOPE\_$PROTEIN-all-IDs.tab
                      head -n 9999 results/$REF/$DEPTH/$EPITOPE/$PROTEIN/ID/*.tab > $CATRAW
                      # perl $CALCULATE $CATRAW $PROTEIN RANDOM_SEQ_$LENGTH >> $SUMMARY/$DEPTH\_$EPITOPE\_$PROTEIN\_summary.txt
                      python $TALLY -i results/$REF/$DEPTH/$EPITOPE/$PROTEIN/ID/ -e $EPITOPE -t $PROTEIN >> $SUMMARY
                done
                #break
        done
        #exit
done
