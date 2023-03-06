#!/bin/bash

# seqkit
# bowtie

module load gcc
module load samtools
module load anaconda3
source activate my-genopipe-env

INSERT=../SyntheticEpitope/scripts/insert_FASTA_into_Genome.pl
RTAG=../db/hg19_EpiID/FASTA_tag/Tag_DB/LAP-tag.fa

GDIR=results/BrowserData/Genomes
[ -d $GDIR ] || mkdir -p $GDIR

# Make two synthetic genomes to align to (opposite strand for rev complement)

# ID3|chr1:23885453-XXX:-
SGENOME=$GDIR/hg19_ID3-Cterm_LAP-tag.fa
perl $INSERT ../input/hg19.fa chr1:23885453:+ $RTAG $SGENOME
bowtie2-build $SGENOME $SGENOME

# NR4A1|chr12:XXXXX-52452725:+
SGENOME=$GDIR/hg19_NR4A1-Cterm_LAP-tag.fa
perl $INSERT ../input/hg19.fa chr12:52452725:- $RTAG $SGENOME
bowtie2-build $SGENOME $SGENOME
