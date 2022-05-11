#!/bin/bash

# Required software:
# python v2.15 with pysam

usage()
{
    echo 'identify-Strain.sh -i /path/to/BAM -g /path/to/genome/fASTA -v /path/to/VCF/files -o /path/to/output [ -s intSeed (default=None) ]'
    echo 'eg: bash identify-Strain.sh -i /input -g /sacCer3.fa -v /sacCer3_VCF -o /output -s 5'
    exit
}

if [ "$#" -ne 8 && "$#" -ne 10 ]; then
    usage
fi

SEED=""

while getopts ":i:g:v:o:s:" IN; do
    case "${IN}" in
        i)
            INPUT=${OPTARG}
            ;;
        g)
            GENOME=${OPTARG}
            ;;
        v)
            VCF=${OPTARG}
            ;;
	o)
            OUTPUT=${OPTARG}
            ;;
	s)
            SEED=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${INPUT}" ] || [ -z "${GENOME}" ] || [ -z "${VCF}" ] || [ -z "${OUTPUT}" ]; then
    usage
fi

echo "Input folder = ${INPUT}"
echo "Genome FASTA file = ${GENOME}"
echo "VCF folder = ${VCF}"
echo "Output folder = ${OUTPUT}"
echo "Seed value = ${SEED}"

LOCAL=$(pwd)
cd $INPUT
for BAM in *.bam
do

	SAMPLE=`basename $BAM`
	echo $SAMPLE

	if [[ $SEED -eq "" ]]; then
		python $LOCAL/strainScripts/detect_strain_BAM.py -b $BAM -g $GENOME -v $VCF -o $OUTPUT/$SAMPLE\_strain.tab
	else
		python $LOCAL/strainScripts/detect_strain_BAM.py -b $BAM -g $GENOME -v $VCF -o $OUTPUT/$SAMPLE\_strain.tab -s $SEED
	fi
done
