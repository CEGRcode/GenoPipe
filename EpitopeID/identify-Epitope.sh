#!/bin/bash

# Required software:
# BWA v0.7.14+
# samtools v1.7+
# bedtools v2.26+
# perl5
# GNU grep (BSD grep on MacOSX is >10X slower)

usage()
{
    echo 'identify-Epitope.sh -i /path/to/FASTQ -o /path/to/output -d /path/to/genome/database -t <Threads - Default 1>'
    echo 'eg: bash identify-Epitope.sh -i /input -o /output -d /sacCer3_EpiID -t 2'
    exit
}

if [ "$#" -ne 6 ] && [ "$#" -ne 8 ]; then
    usage
fi

THREAD=1

while getopts ":i:o:d:t:" IN; do
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
        t)
            THREAD=${OPTARG}
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
echo "Threads = ${THREAD}"

# Check if ALL_TAG.fa bwa index exists, creates if it doesn't
if [ ! -f $DATABASE/tag_FASTA/ALL_TAG.fa.amb ] || [ ! -f $DATABASE/tag_FASTA/ALL_TAG.fa.ann ] || [ ! -f $DATABASE/tag_FASTA/ALL_TAG.fa.bwt ] || [ ! -f $DATABASE/tag_FASTA/ALL_TAG.fa.pac ] || [ ! -f $DATABASE/tag_FASTA/ALL_TAG.fa.sa ]; then
        echo "Building TAG index..."
        bwa index $DATABASE/tag_FASTA/ALL_TAG.fa
fi

# Check if gene.fasta bwa index exists, creates if it doesn't
if [ ! -f $DATABASE/gene_FASTA/gene.fasta.amb ] || [ ! -f $DATABASE/gene_FASTA/gene.fasta.ann ] || [ ! -f $DATABASE/gene_FASTA/gene.fasta.bwt ] || [ ! -f $DATABASE/gene_FASTA/gene.fasta.pac ] || [ ! -f $DATABASE/gene_FASTA/gene.fasta.sa ]; then
        echo "Building ORF index..."
        bwa index $DATABASE/gene_FASTA/gene.fasta
fi

cd $INPUT
for READ1 in *R1*.fastq.gz
do

	SAMPLE=$(echo $READ1 | awk -F"." '{print $1}')
	echo $SAMPLE
	READ2="${READ1/_R1/_R2}"
	
	# Make temporary directory for alignment intermediary files
	mkdir -p $OUTPUT/$SAMPLE
	
	if [ -f $READ2 ]; then
	        echo "Paired-end data detected, aligning to epitope sequence and gene bodies..."
	else
	        echo "Single-end data detected, aligning only to epitope sequence..."
	fi
	
	# Align sequence reads to all tag database
	bwa mem -t $THREAD $DATABASE/tag_FASTA/ALL_TAG.fa $READ1 | samtools view -h -F 4 -S - > $OUTPUT/$SAMPLE/align1.sam
	# Determine which epitope sequences had reads aligned to them
        grep -v "^@" $OUTPUT/$SAMPLE/align1.sam | cut -f1,3 > $OUTPUT/$SAMPLE/epitope-se.out
	
	if [ -f $READ2 ]; then
		bwa mem -t $THREAD $DATABASE/tag_FASTA/ALL_TAG.fa $READ2 | samtools view -h -F 4 -S - > $OUTPUT/$SAMPLE/align2.sam

		# Determine which epitope sequences had reads aligned to them
                grep -v "^@" $OUTPUT/$SAMPLE/align2.sam | cut -f1,3 >> $OUTPUT/$SAMPLE/epitope-se.out

		# Get unique Illumina read ID for each aligned read
		cut -f 1 $OUTPUT/$SAMPLE/align1.sam | grep -v "^@" | uniq > $OUTPUT/$SAMPLE/reads1
		cut -f 1 $OUTPUT/$SAMPLE/align2.sam | grep -v "^@" | uniq > $OUTPUT/$SAMPLE/reads2
		
		# Remove reads pairs that BOTH mapped to the epitope sequence
		grep -Fxv -f $OUTPUT/$SAMPLE/reads1 $OUTPUT/$SAMPLE/reads2 > $OUTPUT/$SAMPLE/uniqueread2
		grep -Fxv -f $OUTPUT/$SAMPLE/reads2 $OUTPUT/$SAMPLE/reads1 > $OUTPUT/$SAMPLE/uniqueread1
		
		# Get the FASTA sequence for the mate-read to the read that mapped to the epitope and remove duplicated paired-end reads	
		gunzip -c $READ1 | grep -h -A 1 -F -f $OUTPUT/$SAMPLE/uniqueread1 - > $OUTPUT/$SAMPLE/read1.fa
		gunzip -c $READ2 | grep -h -A 1 -F -f $OUTPUT/$SAMPLE/uniqueread1 - > $OUTPUT/$SAMPLE/read2.fa
		sed -i -e 's/^@/>/g ; /\--/d' $OUTPUT/$SAMPLE/read1.fa $OUTPUT/$SAMPLE/read2.fa
		perl $DATABASE/epiScripts/uniq_PE_FASTQ.pl $OUTPUT/$SAMPLE/read2.fa $OUTPUT/$SAMPLE/read1.fa > $OUTPUT/$SAMPLE/tag-reads.fa
                gunzip -c $READ2 | grep -h -A 1 -F -f $OUTPUT/$SAMPLE/uniqueread2 - > $OUTPUT/$SAMPLE/read2.fa
                gunzip -c $READ1 | grep -h -A 1 -F -f $OUTPUT/$SAMPLE/uniqueread2 - > $OUTPUT/$SAMPLE/read1.fa
                sed -i -e 's/^@/>/g ; /\--/d' $OUTPUT/$SAMPLE/read2.fa $OUTPUT/$SAMPLE/read1.fa
                perl $DATABASE/epiScripts/uniq_PE_FASTQ.pl $OUTPUT/$SAMPLE/read1.fa $OUTPUT/$SAMPLE/read2.fa >> $OUTPUT/$SAMPLE/tag-reads.fa
	
		# Align TAG-aligned mates to ORF sequences
		bwa mem -t $THREAD $DATABASE/gene_FASTA/gene.fasta $OUTPUT/$SAMPLE/tag-reads.fa | samtools view -h -F 4 -S - > $OUTPUT/$SAMPLE/orf.sam
		# Blacklist filter out artefactual or unwanted gene alignments
		samtools view -Shb $OUTPUT/$SAMPLE/orf.sam > $OUTPUT/$SAMPLE/orf.bam
		bedtools intersect -v -abam $OUTPUT/$SAMPLE/orf.bam -b $DATABASE/blacklist_filter/blacklist.bed | samtools view -h - > $OUTPUT/$SAMPLE/orf.sam
		
		# Determine which epitope sequences had reads aligned to them
		grep -v "^@" $OUTPUT/$SAMPLE/align1.sam | cut -f1,3 > $OUTPUT/$SAMPLE/epitope-pe.out
		grep -v "^@" $OUTPUT/$SAMPLE/align2.sam | cut -f1,3 >> $OUTPUT/$SAMPLE/epitope-pe.out
	fi

        # Parse single-end epitope alignments to final table
        PARSE=$DATABASE/epiScripts/count_raw_epitope.pl
        perl $PARSE $OUTPUT/$SAMPLE/epitope-se.out $OUTPUT/$SAMPLE\-TAG-ID.out

	if [ -f $READ2 ]; then
		# Parse paired-end epitope alignments to final table
		PARSE=$DATABASE/epiScripts/count_tagAlign.pl
		perl $PARSE $OUTPUT/$SAMPLE/epitope-pe.out $OUTPUT/$SAMPLE/orf.sam $OUTPUT/$SAMPLE\-TAG-ID.out
	fi

	# Clean up temporary files
	rm -r $OUTPUT/$SAMPLE

done
