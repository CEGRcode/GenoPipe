
#GENE=REB1
#OUTPUT=sacCer3_Reb1_10M
GENE=$1
OUTPUT=$2

for i in {1..1000}
do
	IDFILE=results/$OUTPUT/ID/Simulation_$i\_deletion.tab
	[ -f $IDFILE ] || continue
	
	NLINE=`wc -l $IDFILE | awk '{print $1}'`
	FIND=`grep $GENE $IDFILE | wc -l | awk '{print $1}'`
	
	if [[ " $NLINE " =~ " 1 " ]] && [[ " $FIND " =~ " 1 " ]];
	then
		echo $NLINE $FIND $IDFILE
	fi
done

