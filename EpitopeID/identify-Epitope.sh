#!/bin/bash

# Required software:
# bowtie2 v2.4.5
# samtools v1.7+
# bedtools v2.26+
# perl5
# python3 with scipy
# GNU grep (BSD grep on MacOSX is >10X slower)

# Unused Software (no-longer required)
# BWA v0.7.14+

usage()
{
	echo 'identify-Epitope.sh -i /path/to/FASTQ -o /path/to/output -d /path/to/genome/database [-t <Threads - Default 1>] [-p <Pvalue - Default 0.05>]'
	echo 'eg: bash identify-Epitope.sh -i /input -o /output -d /sacCer3_EpiID -t 2 -p 0.1'
	exit
}

if [ "$#" -ne 6 ] && [ "$#" -ne 8 ]; then
	usage
fi

PVALUE=0.05
THREAD=1

while getopts ":i:o:d:t:p:" IN; do
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
		p)
			PVALUE=${OPTARG}
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
echo "P-value = ${PVALUE}"

# Check if ALL_TAG.fa exists, exit otherwise
if [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa ]; then
	echo "Epitope FASTA sequences not found: $DATABASE/FASTA_tag/ALL_TAG.fa"
	exit
fi

# Check if GENOME FASTA exists, exit otherwise
if [ ! -f $DATABASE/FASTA_genome/genome.fa ]; then
	echo "Genome FASTA sequence not found: $DATABASE/FASTA_genome/genome.fa"
	exit
fi

# Check if ALL_TAG.fa aligner index exists, creates if it doesn't
#if [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.amb ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.ann ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.bwt ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.pac ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.sa ]; then
if [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.1.bt2 ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.2.bt2 ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.3.bt2 ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.4.bt2 ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.rev.1.bt2 ] || [ ! -f $DATABASE/FASTA_tag/ALL_TAG.fa.rev.2.bt2 ]; then
	echo "Building TAG index..."
	#bwa index $DATABASE/FASTA_tag/ALL_TAG.fa
	bowtie2-build $DATABASE/FASTA_tag/ALL_TAG.fa $DATABASE/FASTA_tag/ALL_TAG.fa
fi

# Check if genome.fa aligner index exists, creates if it doesn't
#if [ ! -f $DATABASE/FASTA_genome/genome.fa.amb ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.ann ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.bwt ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.pac ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.sa ]; then
if [ ! -f $DATABASE/FASTA_genome/genome.fa.1.bt2 ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.2.bt2 ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.3.bt2 ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.4.bt2 ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.rev.1.bt2 ] || [ ! -f $DATABASE/FASTA_genome/genome.fa.rev.2.bt2 ]; then
	echo "Building ORF index..."
	#bwa index $DATABASE/FASTA_genome/genome.fa
	bowtie2-build $DATABASE/FASTA_genome/genome.fa $DATABASE/FASTA_genome/genome.fa
fi

# Folder containing EpitopeID scripts must be located in the same directory as the identify-Epitope.sh shell script
LOCAL=$(pwd)

cd $INPUT
for READ1 in *_R1*.fastq.gz
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
	#bwa mem -t $THREAD $DATABASE/FASTA_tag/ALL_TAG.fa $READ1 | samtools view -h -F 4 -S - > $OUTPUT/$SAMPLE/align1.sam
	bowtie2 --local -p $THREAD -x $DATABASE/FASTA_tag/ALL_TAG.fa -U $READ1 | samtools view -h -F 4 -S - > $OUTPUT/$SAMPLE/align1.sam
	# Determine which epitope sequences had reads aligned to them
	grep -v "^@" $OUTPUT/$SAMPLE/align1.sam | cut -f1,3 > $OUTPUT/$SAMPLE/epitope-se.out

	if [ -f $READ2 ]; then
		#bwa mem -t $THREAD $DATABASE/FASTA_tag/ALL_TAG.fa $READ2 | samtools view -h -F 4 -S - > $OUTPUT/$SAMPLE/align2.sam
		bowtie2 --local -p $THREAD -x $DATABASE/FASTA_tag/ALL_TAG.fa -U $READ2 | samtools view -h -F 4 -S - > $OUTPUT/$SAMPLE/align2.sam

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
		perl $LOCAL/epiScripts/uniq_PE_FASTQ.pl $OUTPUT/$SAMPLE/read2.fa $OUTPUT/$SAMPLE/read1.fa > $OUTPUT/$SAMPLE/tag-reads.fa
		gunzip -c $READ2 | grep -h -A 1 -F -f $OUTPUT/$SAMPLE/uniqueread2 - > $OUTPUT/$SAMPLE/read2.fa
		gunzip -c $READ1 | grep -h -A 1 -F -f $OUTPUT/$SAMPLE/uniqueread2 - > $OUTPUT/$SAMPLE/read1.fa
		sed -i -e 's/^@/>/g ; /\--/d' $OUTPUT/$SAMPLE/read2.fa $OUTPUT/$SAMPLE/read1.fa
		perl $LOCAL/epiScripts/uniq_PE_FASTQ.pl $OUTPUT/$SAMPLE/read1.fa $OUTPUT/$SAMPLE/read2.fa >> $OUTPUT/$SAMPLE/tag-reads.fa

		# Align TAG-aligned mates to genome, requiring mapping quality score of at least 5
		#bwa mem -t $THREAD $DATABASE/FASTA_genome/genome.fa $OUTPUT/$SAMPLE/tag-reads.fastq | samtools view -F 4 -q 5 -Shb - > $OUTPUT/$SAMPLE/orf.bam
		bowtie2 --local -f -p $THREAD -x $DATABASE/FASTA_genome/genome.fa -U $OUTPUT/$SAMPLE/tag-reads.fa | samtools view -F 4 -q 5 -Shb - > $OUTPUT/$SAMPLE/orf.bam

		# Blacklist filter out artefactual or unwanted gene alignments
		bedtools intersect -v -abam $OUTPUT/$SAMPLE/orf.bam -b $DATABASE/blacklist_filter/blacklist.bed > $OUTPUT/$SAMPLE/orf_filter.bam
		# Intersect genomic reads with matched epitope reads with genomic annotation
		bedtools intersect -wb -abam $OUTPUT/$SAMPLE/orf_filter.bam -b $DATABASE/annotation/genome_annotation.gff.gz -bed > $OUTPUT/$SAMPLE/align-pe.out
		# Compress BAM file read length to 1 bp to prevent multi-counting across BIN junctures
		perl $LOCAL/epiScripts/filter_intersect_by_FivePrime.pl $OUTPUT/$SAMPLE/align-pe.out $OUTPUT/$SAMPLE/align-pe_filter.out
	fi

	# Parse single-end epitope alignments to final table
	perl $LOCAL/epiScripts/count_raw_epitope.pl $OUTPUT/$SAMPLE/epitope-se.out $OUTPUT/$SAMPLE/SE_table.out

	if [ -f $READ2 ]; then
		# Parse paired-end epitope alignments to final table
		perl $LOCAL/epiScripts/sum_PE_epitope-alignment.pl $OUTPUT/$SAMPLE/epitope-se.out $OUTPUT/$SAMPLE/align-pe_filter.out $OUTPUT/$SAMPLE/PE_table.out

		echo "Calculating significance..."
		read EPICOUNT ID <<< $(wc -l $OUTPUT/$SAMPLE/epitope-se.out)
		samtools view -H $OUTPUT/$SAMPLE/orf.bam > $OUTPUT/$SAMPLE/sam-header.txt
		GENOMESIZE="$(perl $LOCAL/epiScripts/sum_GenomeSize.pl $OUTPUT/$SAMPLE/sam-header.txt)"

		python $LOCAL/epiScripts/calculate_EpitopeSignificance.py -t $OUTPUT/$SAMPLE/PE_table.out -p $PVALUE -c $EPICOUNT -s $GENOMESIZE -o $OUTPUT/$SAMPLE/PE_sig.out
		cat $OUTPUT/$SAMPLE/SE_table.out $OUTPUT/$SAMPLE/PE_sig.out > $OUTPUT/$SAMPLE\-ID.tab
	else
		cat $OUTPUT/$SAMPLE/SE_table.out > $OUTPUT/$SAMPLE\-ID.tab
	fi

	# Clean up temporary files
	rm -r $OUTPUT/$SAMPLE
done
