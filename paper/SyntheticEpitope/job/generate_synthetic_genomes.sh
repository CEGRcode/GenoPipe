#!/bin/bash

# This script makes genomes to simulate from. Two yeast and two human genomes each with a 500bp Random epitope tag squence at the specificed loci

INSERT=scripts/insert_FASTA_into_Genome.pl
RTAG=../input/RAND_500.fa

[ -d synthetic_genome ] || mkdir synthetic_genome

for RLEN in "500" "100" "50" "20";
do
	RTAG=../input/RAND_$RLEN.fa
	# sacCer3 synthetic genomes
	perl $INSERT ../input/sacCer3.fa chr2:334385:- $RTAG synthetic_genome/sacCer3_Reb1-Cterm_R$RLEN.fa
	perl $INSERT ../input/sacCer3.fa chr14:241688:+ $RTAG synthetic_genome/sacCer3_Rap1-Nterm_R$RLEN.fa

	perl $INSERT ../input/sacCer3.fa chr16:711138:+ $RTAG synthetic_genome/sacCer3_Sua7-Cterm_R$RLEN.fa
	perl $INSERT ../input/sacCer3.fa chr3:205397:- $RTAG synthetic_genome/sacCer3_Taf2-Nterm_R$RLEN.fa

	perl $INSERT ../input/sacCer3.fa chr7:617516:- $RTAG synthetic_genome/sacCer3_Spt4-Cterm_R$RLEN.fa
	perl $INSERT ../input/sacCer3.fa chr2:401253:- $RTAG synthetic_genome/sacCer3_Spt7-Nterm_R$RLEN.fa

	perl $INSERT ../input/sacCer3.fa chr7:998188:+ $RTAG synthetic_genome/sacCer3_Gcn5-Cterm_R$RLEN.fa
	perl $INSERT ../input/sacCer3.fa chr7:368753:+ $RTAG synthetic_genome/sacCer3_Hsf1-Nterm_R$RLEN.fa

	perl $INSERT ../input/sacCer3.fa chr2:586547:- $RTAG synthetic_genome/sacCer3_Fzo1-Cterm_R$RLEN.fa
	perl $INSERT ../input/sacCer3.fa chr16:454990:- $RTAG synthetic_genome/sacCer3_Lge1-Nterm_R$RLEN.fa

	# hg19 synthetic genomes
	perl $INSERT ../input/hg19.fa chr16:67650679:+ $RTAG synthetic_genome/hg19_CTCF-Nterm_R$RLEN.fa
	perl $INSERT ../input/hg19.fa chr3:184079502:+ $RTAG synthetic_genome/hg19_POLR2H-Nterm_R$RLEN.fa

	perl $INSERT ../input/hg19.fa chrX:70338605:+ $RTAG synthetic_genome/hg19_MED12-Nterm_R$RLEN.fa
	perl $INSERT ../input/hg19.fa chr14:100705582:+ $RTAG synthetic_genome/hg19_YY1-Nterm_R$RLEN.fa

	perl $INSERT ../input/hg19.fa chr1:161013065:- $RTAG synthetic_genome/hg19_USF1-Nterm_R$RLEN.fa
	perl $INSERT ../input/hg19.fa chr21:27113910:+ $RTAG synthetic_genome/hg19_GABPA-Nterm_R$RLEN.fa

	perl $INSERT ../input/hg19.fa chr6:152129048:+ $RTAG synthetic_genome/hg19_ESR1-Nterm_R$RLEN.fa
	perl $INSERT ../input/hg19.fa chr14:38064177:- $RTAG synthetic_genome/hg19_FOXA1-Nterm_R$RLEN.fa

	perl $INSERT ../input/hg19.fa chr7:155604816:- $RTAG synthetic_genome/hg19_SHH-Nterm_R$RLEN.fa
	perl $INSERT ../input/hg19.fa chr22:41489009:+ $RTAG synthetic_genome/hg19_EP300-Nterm_R$RLEN.fa
done

for FASTA in `ls synthetic_genome/*.fa`;
do
	samtools faidx $FASTA
done
