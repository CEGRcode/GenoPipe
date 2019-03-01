# This script merges N-number of coverage files generated from calculate_genomic_coverage.sh into a
# single file 'mappability.tab' which is the input for the detect_deletion_BAM.py script

for file in *.out
do
	cut -f1 $file > header
	break
done

for file in *.out
do
	cut -f2 $file > $file.temp
done

paste header *.temp > mappability.tab
rm -f header *.temp
