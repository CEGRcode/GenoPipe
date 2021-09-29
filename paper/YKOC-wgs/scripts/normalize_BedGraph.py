from os import listdir
from os.path import isfile, join
import sys
import re
import argparse

# Python 3.6+
# relies on dict insertion order

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Parse bedgraph file to total tag normalize.')
	
	parser.add_argument('-i','--bedgraph', metavar='bedgraph_fn', dest='bedgraph_fn', required=True, help='the BedGraph file')
	
	args = parser.parse_args()
	return(args)

# chr1	0	1	14
# chr1	1	2	20
# chr1	2	3	25
# chr1	3	5	27
def parse_bedgraph(bg_fn):
	normalize_factor = 0.0
	reader = open(bg_fn, 'r')
	for line in reader:
		if(line.find('#')==0):
			continue
		tokens = line.strip().split('\t')
		start = int(tokens[1])
		stop = int(tokens[2])
		score = float(tokens[3])
		normalize_factor += score*(stop-start)
	reader.close()
	return(normalize_factor)

if __name__ == "__main__":
	'''Normalize total tag of Bedgraph'''
	# Get params
	args = getParams()
	BG_FILE = args.bedgraph_fn
	
	# Pileup BedGraph in CDT interval, then normalize
	if(isfile(BG_FILE)):
		normalize_factor = parse_bedgraph(BG_FILE)
		sys.stderr.write( "Normalization Factor: %i\n" % (normalize_factor))
		reader = open(BG_FILE,'r')
		for line in reader:
			tokens = line.strip().split("\t")
			tokens[3] = str(float(tokens[3])/normalize_factor*1000000)
			sys.stdout.write("\t".join(tokens) + "\n")
		reader.close()
		sys.stderr.write("Normalization Complete for %s\n" % (BG_FILE))
	else:
		sys.stderr.write("File does not exist: %s\n" % (BG_FILE))
