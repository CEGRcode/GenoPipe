#!/bin/bash

# Required software:
# python v2.15 with scipy

usage()
{
    echo 'identify-Deletion.sh -i /path/to/BAM -o /path/to/output -d /path/to/genome/database'
    echo 'eg: bash identify-Deletion.sh -i /input -o /output -d /sacCer3_Del'
    exit
}

if [ "$#" -ne 6 ]; then
    usage
fi

while getopts ":i:o:d:" IN; do
    case "${IN}" in
        i)
            INPUT=${OPTARG}
            ;;
        o)
            OUTPUT=${OPTARG}
            ;;
	d)
	    DATABASE=${OPTARG}
	    ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${INPUT}" ] || [ -z "${OUTPUT}" ] || [ -z "$DATABASE" ]; then
    usage
fi

echo "Input folder = ${INPUT}"
echo "Output folder = ${OUTPUT}"
echo "Databse folder = ${DATABASE}"

LOCAL=$(pwd)
cd $INPUT
for BAM in *.bam
do

	SAMPLE=$(echo $BAM | awk -F"." '{print $1}')
	echo $SAMPLE
	python2 $LOCAL/delScripts/detect_deletion_BAM.py -b $BAM -c $DATABASE/genomic_coord/coord.bed -m $DATABASE/mappability/mappability.tab -o $OUTPUT/$SAMPLE\_deletion.tab -l -2 -M 0.25

done
