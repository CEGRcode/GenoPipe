#!/bin/bash

# sacCer3 VCF database
# This script pulls the variant annotations from SGD in VCF format and parses
# it into a format using arabic numeraled chromosomes. Move the resulting
# database into the desired directory when complete.

#Downloaded from SGD
#Song G., et al., 2014. AGAPE (Automated Genome Analysis PipelinE) for Pan-Genome
#Analysis of Saccharomyces cerevisiae.

# Required software:
# bedtools
# wget
# gzip
# Perl 5.18+

usage()
{
    echo 'generate_sacCer3_VariantDB.sh'
    echo 'eg: bash generate_sacCer3_VariantDB.sh'
    exit
}

[ -d sacCer3_VCF/full_VCF ] || mkdir -p sacCer3_VCF/full_VCF

# Download Strain VCF files
BASEURL=http://sgd-archive.yeastgenome.org/published_datasets/Song_2015_PMID_25781462
wget $BASEURL/BY4741_Stanford_2014_JRIS00000000/BY4741.gatk.vcf.gz
wget $BASEURL/BY4742_Stanford_2014_JRIR00000000/BY4742.gatk.vcf.gz
wget $BASEURL/CEN.PK2-1Ca_Stanford_2014_JRIV01000000/CEN.PK2-1Ca.gatk.vcf.gz
wget $BASEURL/D273-10B_Stanford_2014_JRIY00000000/D273-10B.gatk.vcf.gz
wget $BASEURL/FL100_Stanford_2014_JRIT00000000/FL100.gatk.vcf.gz
wget $BASEURL/JK9-3d_Stanford_2014_JRIZ00000000/JK9-3d.gatk.vcf.gz
wget $BASEURL/RM11-1A_Stanford_2014_JRIP00000000/RM11-1A.gatk.vcf.gz
wget $BASEURL/SEY6210_Stanford_2014_JRIW00000000/SEY6210.gatk.vcf.gz
wget $BASEURL/Sigma1278b-10560-6B_Stanford_2014_JRIQ00000000/Sigma1278b-10560-6B.gatk.vcf.gz
wget $BASEURL/W303_Stanford_2014_JRIU00000000/W303.gatk.vcf.gz
wget $BASEURL/Y55_Stanford_2014_JRIF00000000/Y55.gatk.vcf.gz
mv *.vcf.gz sacCer3_VCF/full_VCF/

# Decompress/unzip all VCF files
for ZVCF in `ls sacCer3_VCF/full_VCF/*.vcf.gz`;
do
	echo "Parsing $ZVCF ..."
	VCF=`basename $ZVCF ".gz"`
	perl parsers/parse_sacCer3_VCF.pl <(gzip -dc $ZVCF) sacCer3_VCF/full_VCF/$VCF
	# clean-up
	rm $ZVCF
done

# Subtract shared variants from sacCer3_VCF
for VCF in `ls sacCer3_VCF/full_VCF/*.vcf`;
do
	STRAIN=`basename $VCF ".gatk.vcf"`
	SUBTRACT=ALL_BUT_$STRAIN

	echo "($STRAIN) Build subtraction file..."
	# write VCF header
	grep '^#' $VCF > $SUBTRACT
	for TVCF in `ls sacCer3_VCF/full_VCF/*.vcf`;
	do
		# append variants that don't match the current strain of interest ($STRAIN)
		if [[ "$TVCF" =~ "$VCF" ]];
		then
			echo "MATCH-$STRAIN"
			continue
		fi
		grep -v '^#' $TVCF >> $SUBTRACT
	done
	echo "Complete"

	echo "($STRAIN) Interesecting shared variants..."
	bedtools intersect -v -header -a $VCF -b $SUBTRACT > sacCer3_VCF/$STRAIN.gatk.vcf
	echo "Complete"

	# clean-up
	rm $SUBTRACT
done

# Check for overlap
#cat *.gatk.vcf |cut -f1-2 |sort |uniq -c |sort -n
