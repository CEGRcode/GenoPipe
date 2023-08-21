#!/bin/bash

# This script makes genomes to simulate from. Two yeast and two human genomes each with a different variant profile based on their respective full sized VCF files

ADDSNPS=scripts/add_VCF_into_Genomic_FASTA.py
SVCF=../db/sacCer3_VCF/full_VCF
HVCF=../db/hg19_VCF
H38VCF=../db/VCF_hg38
[ -d synthetic_genome ] || mkdir synthetic_genome

python $ADDSNPS -f ../input/sacCer3.fa -v $SVCF/RM11-1A.gatk.vcf > synthetic_genome/sacCer3_RM11-1A.fa
python $ADDSNPS -f ../input/sacCer3.fa -v $SVCF/CEN.PK2-1Ca.gatk.vcf > synthetic_genome/sacCer3_CEN.PK2-1Ca.fa

python $ADDSNPS -f ../input/hg19.fa -v $HVCF/K562.vcf > synthetic_genome/hg19_K562.fa
python $ADDSNPS -f ../input/hg19.fa -v $HVCF/HELA.vcf > synthetic_genome/hg19_HELA.fa

python $ADDSNPS -f ../input/hg38.fa -v $H38VCF/K562.vcf > synthetic_genome/hg38_K562.fa
python $ADDSNPS -f ../input/hg38.fa -v $H38VCF/HELA.vcf > synthetic_genome/hg38_HELA.fa

for FASTA in `ls synthetic_genome/*.fa`;
do
	samtools faidx $FASTA
done
