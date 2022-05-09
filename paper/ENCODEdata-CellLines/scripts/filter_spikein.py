import sys,argparse
import pysam

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Remove non-hg19 chromosomes/contigs from the input BAM file.')

	parser.add_argument('-b','--bam', metavar='bamfile', dest='bamfilepath', required=True, help='The BAM data to filter spike in chromosomes from')
	parser.add_argument('-g','--genome', metavar='genome_fasta', dest='genomefilepath', required=True, help='The FASTA with genomic sequence to pull chr names to keep from')
	parser.add_argument('-o','--output', metavar='outfile', dest='outfilepath', required=True, help='The filtered output BAM file path')

	args = parser.parse_args()
	return(args)

if __name__ == '__main__':
	'''Main pileup function'''
	args = getParams()

	# open Genome fasta file
	GENOME = pysam.FastaFile(args.genomefilepath)
	CHR = list(GENOME.references)
	GENOME.close()

	# open BAM files
	samfile = pysam.AlignmentFile(args.bamfilepath, "rb")
	nospike = pysam.AlignmentFile(args.outfilepath,"wb", template=samfile)

	for valid_chr in CHR:
		if(samfile.get_tid(valid_chr) == -1):
			continue
		for read in samfile.fetch(valid_chr):
			nospike.write(read)
	# Close files
	samfile.close()
	nospike.close()
