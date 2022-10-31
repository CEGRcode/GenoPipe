#!/bin/bash

# FIRST CHANGE PATH TO EXECUTE
WRK=/path/to/GenoPipe/paper/SyntheticEpitope
cd $WRK

for EPITOPE in "R500" "R100" "R50" "R20";
do
        for PROTEIN in "Reb1" "Rap1" "Sua7" "Taf2" "Spt4" "Spt7" "Gcn5" "Hsf1" "Fzo1" "Lge1";
        do
                #break
                REF="sacCer3"
                for DEPTH in "10K" "100K" "1M" "10M";
                do
                      head -n 9999 results/$REF/$DEPTH/$EPITOPE/$PROTEIN/ID/*.tab \
                          > results/$REF/$DEPTH/$EPITOPE/$PROTEIN/$REF\_$DEPTH\_$EPITOPE\_$PROTEIN-all-IDs.tab
                done
                break
        done
        exit
done
