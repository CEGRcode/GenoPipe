#!/bin/bash

# Required software:
# BWA v0.7.14+
# samtools v1.7+
# bedtools v2.26+
# perl5
# python v2.15 with scipy
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
if [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.amb ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.ann ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.bwt ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.pac ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.sa ]; then
        echo "Building TAG index..."
        bwa index $DATABASE/FASTA_tag/ALL_TAG.fa
fi

# Check if genome.fa bwa index exists, creates if it doesn't
if [ ! -f $DATABASE/FASTA_genome/genome.fa.amb ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.ann ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.bwt ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.pac ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.sa ]; then
        echo "Building ORF index..."
        bwa index $DATABASE/FASTA_genome/genome.fa
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
	bwa mem -t $THREAD $DATABASE/FASTA_tag/ALL_TAG.fa $READ1 | samtools view -h -F 4 -S - > $OUTPUT/$SAMPLE/align1.sam
	# Determine which epitope sequences had reads aligned to them
        grep -v "^@" $OUTPUT/$SAMPLE/align1.sam | cut -f1,3 > $OUTPUT/$SAMPLE/epitope-se.out
	
	if [ -f $READ2 ]; then
		bwa mem -t $THREAD $DATABASE/FASTA_tag/ALL_TAG.fa $READ2 | samtools view -h -F 4 -S - > $OUTPUT/$SAMPLE/align2.sam

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
	
		# Align TAG-aligned mates to genome, requiring mapping quality score of at least 5
		bwa mem -t $THREAD $DATABASE/FASTA_genome/genome.fa $OUTPUT/$SAMPLE/tag-reads.fa | samtools view -F 4 -q 5 -Shb - > $OUTPUT/$SAMPLE/orf.bam

		# Blacklist filter out artefactual or unwanted gene alignments
		bedtools intersect -v -abam $OUTPUT/$SAMPLE/orf.bam -b $DATABASE/blacklist_filter/blacklist.bed | samtools view -Shb - > $OUTPUT/$SAMPLE/orf_filter.bam
		# Intersect genomic reads with matched epitope reads with genomic annotation
                bedtools intersect -wb -abam $OUTPUT/$SAMPLE/orf_filter.bam -b $DATABASE/annotation/genome_annotation.gff.gz -bed > $OUTPUT/$SAMPLE/align-pe.out
                # Compress BAM file read length to 1 bp to prevent multi-counting across BIN junctures
                perl $DATABASE/epiScripts/filter_intersect_by_FivePrime.pl $OUTPUT/$SAMPLE/align-pe.out $OUTPUT/$SAMPLE/align-pe_filter.out

	fi

        # Parse single-end epitope alignments to final table
        perl $DATABASE/epiScripts/count_raw_epitope.pl $OUTPUT/$SAMPLE/epitope-se.out $OUTPUT/$SAMPLE/SE_table.out

	if [ -f $READ2 ]; then
		# Parse paired-end epitope alignments to final table
		perl $DATABASE/epiScripts/sum_PE_epitope-alignment.pl $OUTPUT/$SAMPLE/epitope-se.out $OUTPUT/$SAMPLE/align-pe_filter.out $OUTPUT/$SAMPLE/PE_table.out

	        echo "Calculating significance..."
	        read EPICOUNT ID <<< $(wc -l $OUTPUT/$SAMPLE/epitope-se.out)
	        samtools view -H $OUTPUT/$SAMPLE/orf.bam > sam-header.txt
	        GENOMESIZE="$(perl $DATABASE/epiScripts/sum_GenomeSize.pl sam-header.txt)"

		PVALUE=0.05
		python2 $DATABASE/epiScripts/calculate_EpitopeSignificance.py -t $OUTPUT/$SAMPLE/PE_table.out -p $PVALUE -c $EPICOUNT -s $GENOMESIZE -o $OUTPUT/$SAMPLE/PE_sig.out
		cat $OUTPUT/$SAMPLE/SE_table.out $OUTPUT/$SAMPLE/PE_sig.out > $OUTPUT/$SAMPLE\-ID.tab
	else
		cat $OUTPUT/$SAMPLE/SE_table.out > $OUTPUT/$SAMPLE\-ID.tab
	fi

	# Clean up temporary files
	rm -r $OUTPUT/$SAMPLE

done
