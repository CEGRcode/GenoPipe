import sys, getopt, csv, os
import pysam

usage = """
Usage:
This script calculates the coverage of the supplied intervals given a BAM file
aligned with uniformly tiling reads

Example: python3 calculate_region_coverage.py -f NGS.bam -c Interval.bed
"""

READLENGTH = 0

def iterateBAM(bam, bed):
	# open BAM file
	bamfile = pysam.AlignmentFile(bam, "rb")
	# open BED file
	file = open(bed, "r")
	
	# Array to hold genomic ID and mappability score
	SCORE = []
	
	# Keep track of read length from BAM file
	global READLENGTH
	
	# Iterate bed coord file, getting tag counts across interval
	for line in file:
		bedline = line.rstrip().split("\t")
		intervalID = bedline[3]
		intervalSize = (int(bedline[2]) - int(bedline[1])) + 1
		if intervalSize < 1:
			try:
				sys.exit(-1)
			except SystemExit:
				print(line.rstrip()+"|Interval coordinates invalid\n")
		else:
			try:
				intervalCount = 0
				for read in bamfile.fetch(bedline[0], int(bedline[1]), int(bedline[2])):
					intervalCount = intervalCount + 1
					if read.query_length > READLENGTH:
						READLENGTH = read.query_length
				MAP = float(intervalCount) / float(intervalSize)
#				print intervalID + "\t" + str(intervalSize) + "\t" + str(intervalCount)
				SCORE.append((intervalID,MAP))
			except (IndexError, ValueError):
				try:
					sys.exit(-1)
				except SystemExit:
					print(line.rstrip()+"|Not present in BAM File\n")
	# Close files
	file.close()
	bamfile.close()
	
	return SCORE

def validateBAM(bam):
	samfile = pysam.AlignmentFile(bam, "rb")
	index = samfile.check_index()
	samfile.close()
	return index

# Main program which takes in input parameters
if __name__ == '__main__':
	if len(sys.argv) < 2 or not sys.argv[1].startswith("-"): sys.exit(usage)
	
	#python2 $COVERAGE -b ALIGN_GENOME.bam -c $COORD -o $READLENGTH\bp_Cov.out
        # get arguments
	optlist, alist = getopt.getopt(sys.argv[1:], 'hb:c:r:o:')
	for opt in optlist:
		if opt[0] == "-h": sys.exit(usage)
		elif opt[0] == "-b": BAM = opt[1]
		elif opt[0] == "-c": BED = opt[1]
		elif opt[0] == "-r": RES = float(opt[1])
		elif opt[0] == "-o": OUT = opt[1]
		else: sys.exit(usage)
	
	# Validate BAM file
	if(not validateBAM(BAM)):
		try:
			sys.exit(-1)
		except SystemExit:
			print("ERROR!!!\tNo BAM index detected.\n")
	
	print("BAM file: %s" % BAM)
	print("BED file: %s" % BED)
	print("Resolution: %f" % RES)
	print("Output file: %s" % OUT)
	
	# Load BED file and calculate coverage using BAM file
	SCORE = iterateBAM(BAM, BED)
	
	output = open(OUT, "w")
	output.write("GeneID\t" + str(READLENGTH) + "\n")
	for id, score in SCORE:
		score *= RES
		output.write(id + "\t" + str(score) + "\n")
	output.close()
	
	print("Program complete")
