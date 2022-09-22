#!/bin/bash

# FIRST CHANGE PATH TO EXECUTE
WRK=/path/to/GenoPipe/paper/SyntheticEpitope
cd $WRK

TEMPLATE=$WRK/job/depth_template.pbs

N_EXPERIMENTS=`wc -l depth_simulations.txt |awk '{print $1}'`

for (( E=1; E<=$N_EXPERIMENTS; E++ ))
do
	INFO=`sed "${E}q;d" depth_simulations.txt`
	REF=`awk '{print $1}' <(echo $INFO)`
	PROTEIN=`awk '{print $2}' <(echo $INFO)`
	LOCUS=`awk '{print $3}'  <(echo $INFO)`
	EPITOPE=`awk '{print $4}'  <(echo $INFO)`
	TIME="00:10:00"
	[[ "$REF" =~ "hg19" ]] && TIME="02:00:00"

	EXPERIMENTNAME=$PROTEIN-$EPITOPE
	FILENAME=job/run_depth_$E\_$PROTEIN-$LOCUS\_$EPITOPE.pbs
	echo $FILENAME
	sed "s/EXPERIMENTID/${E}/g" $TEMPLATE \
		| sed "s/EXPERIMENTNAME/${EXPERIMENTNAME}/g" \
		| sed "s/TIME/${TIME}/g" \
		> $FILENAME
done
