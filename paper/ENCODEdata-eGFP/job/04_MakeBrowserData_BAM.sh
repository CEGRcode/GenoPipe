
# Before running this script, run EpitopeID with "Clean-up" removal of intermediate files commented out

# Create directory for FASTQ subset results
BDIR=results/BrowserData/BAM
FDIR=results/BrowserData/FASTQ
GDIR=results/BrowserData/Genomes
[ -d $BDIR ] || mkdir $BDIR
[ -d $FDIR ] || mkdir $FDIR

# Align reads for each ENCODE sample with ID3-eGFP and NR4A1-eGFP genomes
for ENCFF in "ENCFF548RTA" "ENCFF671VDI";
do
	FQ=results/FASTQ/$ENCFF
	SFQ=$FDIR/$ENCFF

	# Subset based on read IDs from intermediate files of EpitopeID
	seqkit grep -f <(cat results/ID/$ENCFF\_R1/reads*) $FQ\_R1.fastq.gz > $SFQ\_R1.fastq
	seqkit grep -f <(cat results/ID/$ENCFF\_R1/reads*) $FQ\_R1.fastq.gz > $SFQ\_R2.fastq

	# Align to ID3 and index
	SGENOME=$GDIR/hg19_ID3-Cterm_LAP-tag.fa
	BAM=$BDIR/ID3-Nterm-LAP_$ENCFF
	bowtie2 -x $SGENOME -1 $SFQ\_R1.fastq -2 $SFQ\_R2.fastq | samtools sort -o $BAM.bam
	samtools index $BAM.bam

	# Align to NR4A1 and index
	SGENOME=$GDIR/hg19_NR4A1-Cterm_LAP-tag.fa
	BAM=$BDIR/NR4A1-Nterm-LAP_$ENCFF
	bowtie2 -x $SGENOME -1 $SFQ\_R1.fastq -2 $SFQ\_R2.fastq | samtools sort -o $BAM.bam
	samtools index $BAM.bam
done
