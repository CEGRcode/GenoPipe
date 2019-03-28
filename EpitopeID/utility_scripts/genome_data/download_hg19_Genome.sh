#!/bin/bash
  
# hg19
# This script downloads the hg19 genome from UCSC and indexes it for use with BWA.
# Move the resulting files into the appropriate /pwd/genome/ when complete

# Required software:
# wget
# bwa v0.7.14+
# 
# Optional software:
# twoBitToFa

# Download hg19.2bit
wget -N http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

# Check for the existence of twoBitToFa, download if not present
if ! command -v twoBitToFa; then
	unameOUT="$(uname -s)"

	if [ $unameOUT == "Darwin" ]; then
		wget -N http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa
	else
		wget -N http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
	fi
fi

# Make twoBitToFa globally executable
chmod 777 twoBitToFa
echo "Converting 2bit to fa..."
./twoBitToFa hg19.2bit genome.fa

echo "BWA Indexing genome..."
bwa index genome.fa
echo "Complete"

# Clean up
rm twoBitToFa hg19.2bit
