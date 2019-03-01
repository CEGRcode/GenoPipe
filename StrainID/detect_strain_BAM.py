import sys, getopt, csv, os
import math, random
import pysam

usage = """
Usage:
This script examines an aligned BAM file against a folder containing any number of variant (VCF) files.
The BAM file is checked for the presence of the SNPs in order to determine the strain of origin.

Example: python2 detect_strain_BAM.py -b NGS.bam -g Genome.fasta -v /path/to/VCF/files -o NGS_strain.tab
"""

def buildBackground(bam, gen, count):
        # open BAM file
        samfile = pysam.AlignmentFile(bam, "rb")
	# open Genome fasta file
	GENOME = pysam.FastaFile(gen)

	# Select random reads in BAM file to build SNP background model
        RANDREAD = [random.randint(0,samfile.mapped - 1) for _ in range(count)]
        RANDREAD.sort()

	BACK_REF = BACK_ALT = 0

	it = samfile.fetch()
	# Iterate across the randomly selected read indexes
	for INDEX in RANDREAD:
		# Enumerate the samfile iterator and jump to the index directly
		for index, read in  enumerate(it, INDEX):
			LOC = random.randint(1, read.query_length)
			REFSEQ = GENOME.fetch(read.reference_name, read.reference_start + LOC - 1, read.reference_start + LOC).upper()
			if REFSEQ == read.query_sequence[LOC-1:LOC]:
				BACK_REF = BACK_REF + 1
#			elif random.randint(0,2) == 0: # If not reference sequence, then there is a 1/3rd chance that it matches 'SNP' allele
			else:
				BACK_ALT = BACK_ALT + 1

			# Break out the iterator so we can move to the next index directly with the next loop pass
			break
	# Close files
	GENOME.close()
	samfile.close()

        print "Backround nucleotides examined: " + str(count)
        print "Theoretical SNPs detected: " + str(BACK_ALT)
        print "Reference nucleotides detected: " + str(BACK_REF)

#	sys.exit(1)

	return BACK_ALT, BACK_REF

def iterateVCF(bam, vcf):
	# open BAM file
	samfile = pysam.AlignmentFile(bam, "rb")
       	# open VCF file
	vcffile = pysam.VariantFile(vcf, "r")

	TOTALSNPS = 0
	TOTAL_ALT= 0
	TOTAL_REF = 0

	# Iterate bed coord file, getting tag counts across interval
	for VCFrecord in vcffile.fetch():
		TOTALSNPS = TOTALSNPS + 1
                REFCOUNT = ALTCOUNT = 0

#		print VCFrecord.contig + "\t" + str(VCFrecord.pos) + "\tREF: " + VCFrecord.ref + "\tALT: " + str(VCFrecord.alts)
		for read in samfile.fetch(VCFrecord.contig, VCFrecord.pos, VCFrecord.pos+1):
			index = VCFrecord.pos - read.reference_start
                        if read.query_sequence[index-1:index] in VCFrecord.ref:
                                REFCOUNT = REFCOUNT + 1
                        elif read.query_sequence[index-1:index] in VCFrecord.alts:
                                ALTCOUNT = ALTCOUNT + 1
		TOTAL_ALT = TOTAL_ALT + ALTCOUNT
		TOTAL_REF = TOTAL_REF + REFCOUNT
	# Close files
	samfile.close()
	vcffile.close()

	print "SNPs examined: " + str(TOTALSNPS)
	print "SNPs detected: " + str(TOTAL_ALT)
	print "Reference detected: " + str(TOTAL_REF)
	return TOTAL_ALT,TOTAL_REF

def validateBAM(bam):
	try:
	        samfile = pysam.AlignmentFile(bam, "rb")
		index = samfile.check_index()
		samfile.close()
		return index
	except:
		print "BAM index not detected.\nAttempting to index now...\n"
		pysam.index(str(bam))
	        if not os.path.isfile(bam + ".bai"):
        	        raise RuntimeError("BAM indexing failed, please check if BAM file is sorted")     
                	return False
		print "BAM index successfully generated.\n"
		return True

# Main program which takes in input parameters
if __name__ == '__main__':
        if len(sys.argv) < 2 or not sys.argv[1].startswith("-"): sys.exit(usage)

        BAM = GENOME = VCF = OUT = ""

	# get arguments
        optlist, alist = getopt.getopt(sys.argv[1:], 'hb:g:v:o:')
        for opt in optlist:
                if opt[0] == "-h": sys.exit(usage)
		elif opt[0] == "-b": BAM = opt[1]
		elif opt[0] == "-g": GENOME = opt[1]
                elif opt[0] == "-v": VCF = opt[1]
		elif opt[0] == "-o": OUT = opt[1]
		else: sys.exit(usage)

	if BAM == "":
		print "No BAM file detected!!!"
		sys.exit(usage)
	elif GENOME == "":
		print "No Genomic fasta file detected!!!"
		sys.exit(usage)
	elif not os.path.isdir(VCF):
		print "Invalid path to VCF files!!!"
		sys.exit(usage)

	if OUT == "":
                OUT = os.path.splitext(os.path.basename(BAM))[0] + ".tab"

	# Validate BAM file
	if(not validateBAM(BAM)):
		sys.exit(-1)

	print "BAM file: ",BAM
	print "Genome fasta: ",GENOME
        print "VCF folder: ",VCF
	print "Output file: ",OUT

	# Open output file for writing
        output = open(OUT, "w")
	output.write("\t" + os.path.basename(BAM) + "\n")

	for filename in os.listdir(VCF):
	    if filename.endswith(".vcf"):
		print("\nProcessing: " + filename)
		ALT, REF = iterateVCF(BAM, os.path.join(VCF, filename))
		if(ALT > 0):
			BALT, BREF = buildBackground(BAM, GENOME, (ALT + REF))
			if BALT > 0:
				RATIO = float(ALT) / float(ALT + REF)
				BACK = float(BALT) / float(BALT + BREF)
				LOG2 = math.log(RATIO / BACK) / math.log(2)
				print "Strain log2 enrichment: " + str(LOG2)
				output.write(filename + "\t" + str(LOG2) + "\n")
			else:
				print "Strain log2 enrichment: Inf"
				output.write(filename + "\tInf\n")
		else:
			print "Strain log2 enrichment: NaN"
			output.write(filename + "\tNaN\n")
	output.close()
