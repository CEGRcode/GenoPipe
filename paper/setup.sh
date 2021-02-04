#!/bin/bash

# This script sets up input files shared by synthetic and data studies.
#   -generate random 500bp epitope tag used by SyntheticEpitope
#   -downloads the sacCer3 genome from SGD and indexes it for use with BWA.
#   -downloads the hg19 genome from UCSC and indexes it for use with BWA.

# Required software:
# wget
# Perl 5.18+
# bwa v0.7.14+
# 
# Optional software:
# twoBitToFa

WRK=`pwd`

# Generate random FASTA tag
echo "Create RandomSequence (RTAG)"
RANDTAG=scripts/generate_random_FASTA_sequence.pl
RTAG=input/RANDOM_SEQ.fa
[ -f $RTAG ] || perl $RANDTAG 500 0 $RTAG

# Download Yeast Genome (sacCer3)
YGENOME=input/sacCer3.fa
wget -N https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz
tar -xvzf S288C_reference_genome_R64-1-1_20110203.tgz
echo "Parsing genome..."
YPARSER=../EpitopeID/utility_scripts/genome_data/parsers/parse_sacCer3_Genome_FASTA.pl
perl $YPARSER S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa $YGENOME
echo "Complete"
echo "BWA Indexing genome..."
bwa index $YGENOME
echo "Complete"
rm S288C_reference_genome_R64-1-1_20110203.tgz
rm -r S288C_reference_genome_R64-1-1_20110203/

# Download Human Genome (hg19)
HGENOME=input/hg19.fa
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
./twoBitToFa hg19.2bit $HGENOME
echo "BWA Indexing genome..."
bwa index $HGENOME
echo "Complete"
rm twoBitToFa hg19.2bit
