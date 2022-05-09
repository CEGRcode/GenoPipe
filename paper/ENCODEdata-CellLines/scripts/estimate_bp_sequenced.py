import sys, os
import math, random
import argparse
import pysam

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Sample reads from a BAM file to estimate total bp sequenced')

	parser.add_argument('-b','--bam', metavar='bamfile', dest='bamfilepath', required=True, help='The BAM file to sample and estimate bp sequenced from')

	args = parser.parse_args()
	return(args)

def validateBAM(bam):
	try:
		samfile = pysam.AlignmentFile(bam, "rb")
		index = samfile.check_index()
		samfile.close()
		return index
	except:
		print("BAM index not detected.\nAttempting to index now...\n")
		pysam.index(str(bam))
		if not os.path.isfile(bam + ".bai"):
			raise RuntimeError("BAM indexing failed, please check if BAM file is sorted")
			return False
		print("BAM index successfully generated.\n")
		return True

# Main program which takes in input parameters
if __name__ == '__main__':
	'''Main sampling function'''
	args = getParams()

	RLEN = []
	# Validate BAM file
	if(not validateBAM(args.bamfilepath)):
		sys.exit(-1)
	# Open BAM file
	samfile = pysam.AlignmentFile(args.bamfilepath, "rb")
	# Get num mapped reads
	MAPPED = samfile.mapped

	# Select random reads in BAM file to build SNP background model
	# Sample the larger of 5% of mapped or 1 million reads
	count = max(int((float(MAPPED) * 0.05)), 1000000)
	RANDREAD = [random.randint(0,MAPPED - 1) for _ in range(count)]
	RANDREAD.sort()
	it = samfile.fetch()
	# Iterate across the randomly selected read indexes
	for INDEX in RANDREAD:
		# Enumerate the samfile iterator and jump to the index directly
		for index, read in  enumerate(it, INDEX):
			# Save query/read length
			RLEN.append(read.query_length)
			# Break out the iterator so we can move to the next index directly with the next loop pass
			break
	# Close files
	samfile.close()

	RLEN.sort()
	MEDIAN = RLEN[int(len(RLEN)/2)]
	MEAN = sum(RLEN)/len(RLEN)
	ESTBP = MEDIAN * MAPPED

	sys.stdout.write("\t".join([str(s) for s in [os.path.basename(args.bamfilepath),ESTBP,MAPPED,MEDIAN,max(RLEN),min(RLEN),MEAN]]) + "\n")
