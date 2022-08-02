#!/bin/bash

# This script makes genomes to simulate from. Two yeast and two human genomes each with a 500bp Random epitope tag squence at the specificed loci

INSERT=scripts/insert_FASTA_into_Genome.pl
RTAG=../input/RAND_500.fa

[ -d synthetic_genome ] || mkdir synthetic_genome

for RLEN in "500" "100" "50" "20";
do
	RTAG=../input/RAND_$RLEN.fa
	# sacCer3 synthetic genomes
	perl $INSERT ../input/sacCer3.fa chr14:241688:+ $RTAG synthetic_genome/sacCer3_Rap1-Nterm_R$RLEN.fa
	perl $INSERT ../input/sacCer3.fa chr2:334385:- $RTAG synthetic_genome/sacCer3_Reb1-Cterm_R$RLEN.fa
	# hg19 synthetic genomes
	perl $INSERT ../input/hg19.fa chr16:67650679:+ $RTAG synthetic_genome/hg19_CTCF-Nterm_R$RLEN.fa
	perl $INSERT ../input/hg19.fa chr3:184079502:+ $RTAG synthetic_genome/hg19_POLR2H-Nterm_R$RLEN.fa
done

for FASTA in `ls synthetic_genome/*.fa`;
do
	samtools faidx $FASTA
done
