import sys
import argparse

valid_chr = ["chr1", "chr2",
	"chr3", "chr4",
	"chr5", "chr6",
	"chr7", "chr8",
	"chr9", "chr10",
	"chr11", "chr12",
	"chr13", "chr14",
	"chr15", "chr16",
	"chr17", "chr18",
	"chr19", "chr20",
	"chr21", "chr22",
	"chrM",
	"chrX", "chrY"]

def get_params():
	parser = argparse.ArgumentParser(description="parse hg19 to strip out unwanted chromosomes (haplotype chr) and write new FSTA to STDOUT")
	parser.add_argument(dest="genome_fn", help="input file with two matrices", metavar="GENOME.fa")
	return(parser.parse_args())

if __name__=='__main__':
	args = get_params()
	valid =False
	reader = open(args.genome_fn,'r')
	line = None
	while( line!="" ):
		line = reader.readline().strip()
		if( line[1:] in valid_chr ):
			valid=True
		elif( line.find('>')==0 ):
			valid=False
		if(valid):
			sys.stdout.write(line+'\n')
	reader.close()
