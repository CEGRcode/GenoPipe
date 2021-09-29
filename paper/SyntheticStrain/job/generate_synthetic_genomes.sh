#!/bin/bash

# This script makes genomes to simulate from. Two yeast and two human genomes each with variants from VCF incorporated into sequence

ADDSNPS=scripts/add_VCF_into_Genomic_FASTA.py
SVCF=../db/sacCer3_VCF
HVCF=../db/hg19_VCF
[ -d synthetic_genome ] || mkdir synthetic_genome

python $ADDSNPS -f ../input/sacCer3.fa -v $SVCF/RM11-1A.gatk.vcf > synthetic_genome/sacCer3_RM11-1A.fa
python $ADDSNPS -f ../input/sacCer3.fa -v $SVCF/CEN.PK2-1Ca.gatk.vcf > synthetic_genome/sacCer3_CEN.PK2-1Ca.fa

python $ADDSNPS -f ../input/hg19.fa -v $HVCF/K562.vcf > synthetic_genome/hg19_K562.fa
python $ADDSNPS -f ../input/hg19.fa -v $HVCF/HELA.vcf > synthetic_genome/hg19_HELA.fa
