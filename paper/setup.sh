#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=16gb
#PBS -l walltime=03:00:00
#PBS -A open
#PBS -o run.setup.out
#PBS -e run.setup.err

# This script sets up input files shared by synthetic and data studies.
#   -generate random 500bp epitope tag used by SyntheticEpitope
#   -downloads the sacCer3 genome from SGD and indexes it for use with BWA.
#   -downloads the hg19 genome from UCSC and indexes it for use with BWA.

# Required software:
# wget
# Python 3
# Perl 5.18+
# bwa v0.7.14+
# bowtie v1.2.3
#
# Optional software:
# twoBitToFa

WRK=/path/to/GenoPipe/paper
cd $WRK

module load gcc
module load bwa
module load samtools
module load anaconda3
source activate my-genopipe-env

[ -d $WRK/input ] || mkdir $WRK/input
[ -d $WRK/db ] || mkdir $WRK/db

## Generate random FASTA tag
#echo "Create RandomSequence (RTAG)"
#RANDTAG=scripts/generate_random_FASTA_sequence.pl
#RTAG=input/RANDOM_SEQ.fa
#[ -f $RTAG ] || perl $RANDTAG 500 0 $RTAG

# Download Yeast Genome (sacCer3)
YGENOME=input/sacCer3.fa
if [[ ! -f $YGENOME ]]; then
	echo "**Yeast genome not found, downloading to $YGENOME..."
	wget -N https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz
	tar -xvzf S288C_reference_genome_R64-1-1_20110203.tgz
	echo "Parsing genome..."
	YPARSER=../EpitopeID/utility_scripts/genome_data/parsers/parse_sacCer3_Genome_FASTA.pl
	perl $YPARSER S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa $YGENOME
	echo "BWA Indexing genome..."
	bwa index $YGENOME
	echo "Bowtie2 Indexing genome..."
	bowtie2-build $YGENOME $YGENOME
	rm S288C_reference_genome_R64-1-1_20110203.tgz
	rm -r S288C_reference_genome_R64-1-1_20110203/
fi

# Download Human Genome (hg19)
HGENOME=input/hg19.fa
if [[ ! -f $HGENOME ]]; then
	echo "**Human genome not found, downloading to $HGENOME..."
	wget -N http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
	# Check for the existence of twoBitToFa, download if not present and make globally executable
	if ! command -v twoBitToFa; then
		unameOUT="$(uname -s)"
		if [ $unameOUT == "Darwin" ]; then
			wget -N http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa
		else
			wget -N http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/twoBitToFa
		fi
	fi
	chmod 777 twoBitToFa
	echo "Converting 2bit to fa..."
	./twoBitToFa hg19.2bit $HGENOME.raw
	echo "Strip haplotypes..."
	python scripts/parse_hg19_Genome_FASTA.py $HGENOME.raw > $HGENOME
	echo "BWA Indexing genome..."
	bwa index $HGENOME
	echo "Bowtie2 Indexing genome..."
	bowtie2-build $HGENOME $HGENOME
	echo "Complete"
	rm twoBitToFa hg19.2bit $HGENOME.raw
fi


# Download Human Genome (hg38)
HGENOME=input/hg38.fa
if [[ ! -f $HGENOME ]]; then
	echo "**Human genome not found, downloading to $HGENOME..."
	wget -N http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
	# Check for the existence of twoBitToFa, download if not present and make globally executable
	if ! command -v twoBitToFa; then
		unameOUT="$(uname -s)"
		if [ $unameOUT == "Darwin" ]; then
			wget -N http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa
		else
			wget -N http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/twoBitToFa
		fi
	fi
	chmod 777 twoBitToFa
	echo "Converting 2bit to fa..."
	./twoBitToFa hg38.2bit $HGENOME
	echo "BWA Indexing genome..."
	bwa index $HGENOME
	rm twoBitToFa hg38.2bit
fi

# Build Yeast EpiID with R500
YEPIDB=db/sacCer3_EpiID
[ -d $YEPIDB ] || cp -r ../EpitopeID/sacCer3_EpiID/ $YEPIDB
if [ ! -f $YEPIDB/FASTA_genome/genome.fa ]; then
	echo "**Setup Yeast EpiID genome..."
	cp input/sacCer3.fa $YEPIDB/FASTA_genome/genome.fa
	echo "BWA Indexing genome..."
	bwa index $YEPIDB/FASTA_genome/genome.fa
	echo "Bowtie2 Indexing genome..."
	bowtie2-build $YEPIDB/FASTA_genome/genome.fa $YEPIDB/FASTA_genome/genome.fa
fi
echo "Setup Random Epitope for Yeast EpiID..."
cp input/RAND_500.fa $YEPIDB/FASTA_tag/Tag_DB/
cp input/RAND_100.fa $YEPIDB/FASTA_tag/Tag_DB/
cp input/RAND_50.fa $YEPIDB/FASTA_tag/Tag_DB/
cp input/RAND_20.fa $YEPIDB/FASTA_tag/Tag_DB/
rm $YEPIDB/FASTA_tag/ALL_TAG.fa*
cd $YEPIDB/FASTA_tag/Tag_DB
cat *.fa *.fna *.ffn *.fasta > ALL_TAG.fa
bwa index ALL_TAG.fa
mv ALL_TAG.fa* ../
echo "Complete"
cd $WRK

# Build Human EpiID with R500
HEPIDB=db/hg19_EpiID
[ -d $HEPIDB ] || cp -r ../EpitopeID/hg19_EpiID/ $HEPIDB
if [ ! -f $HEPIDB/FASTA_genome/genome.fa ]; then
	echo "**Setup Human EpiID genome..."
	cp input/hg19.fa $HEPIDB/FASTA_genome/genome.fa
	# echo "BWA Indexing genome..."
	# bwa index $HEPIDB/FASTA_genome/genome.fa
	# echo "Complete"
	echo "Bowtie2 Indexing genome..."
	bowtie2-build $HEPIDB/FASTA_genome/genome.fa $HEPIDB/FASTA_genome/genome.fa
fi
echo "Setup Random Epitope for Human EpiID..."
cp input/RAND_500.fa $HEPIDB/FASTA_tag/Tag_DB/
rm $HEPIDB/FASTA_tag/ALL_TAG.fa*
cd $HEPIDB/FASTA_tag/Tag_DB
cat *.fa *.fna *.ffn *.fasta > ALL_TAG.fa
bwa index ALL_TAG.fa
mv ALL_TAG.fa* ../
echo "Complete"
cd $WRK

# Build Human EpiID with HIV genome
HHIVDB=db/hg19-HIV_EpiID
[ -d $HHIVDB ] || cp -r ../EpitopeID/hg19_EpiID/ $HHIVDB
if [ ! -f $HHIVDB/FASTA_genome/genome.fa ]; then
	echo "**Setup Human EpiID genome..."
	cp input/hg19.fa $HHIVDB/FASTA_genome/genome.fa
	# echo "BWA Indexing genome..."
	# bwa index $HHIVDB/FASTA_genome/genome.fa
	# echo "Complete"
	echo "Bowtie2 Indexing genome..."
	bowtie2-build $HHIVDB/FASTA_genome/genome.fa $HHIVDB/FASTA_genome/genome.fa
	echo "Complete"
fi
echo "Setup HIV genome for Human EpiID..."
cp input/AF324493.2_HIV-1_vector_pNL4-3.fa $HHIVDB/FASTA_tag/Tag_DB/
rm $HHIVDB/FASTA_tag/ALL_TAG.fa*
cd $HHIVDB/FASTA_tag/Tag_DB
cat *.fa *.fna *.ffn *.fasta > ALL_TAG.fa
# bwa index ALL_TAG.fa
bowtie2-build ALL_TAG.fa ALL_TAG.fa
mv ALL_TAG.fa* ../
echo "Complete"

cd $WRK

# Add Yeast Deletion DB to paper/db
cd $WRK/db
ln -s ../../DeletionID/sacCer3_Del
cd $WRK

# Add Yeast & Human StrainID DB to paper/db
cd $WRK/../StrainID
tar -xvf hg38_VCF
tar -xvf hg38_ENCODE
cd utility_scripts
bash generate_hg38_VariantDB.#!/bin/sh
mv hg38_DepMap ../
cd $WRK/db
ln -s ../../StrainID/sacCer3_VCF
ln -s ../../StrainID/hg19_VCF
ln -s ../../StrainID/hg38_VCF
ln -s ../../StrainID/hg38_ENCODE
ln -s ../../StrainID/hg38_DepMap
cd $WRK

# Setup color-space index for yeast genome
# (used by BY4742-chipseq)
bowtie-build -C input/sacCer3.fa input/sacCer3_index
