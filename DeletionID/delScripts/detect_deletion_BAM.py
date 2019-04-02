import sys, getopt, csv, os
import numpy as np
import pysam

usage = """
Usage:
This script determine if any of the supplied intervals are depleted or missing from the BAM file.
Required parameters:
BAM File
BED File
Mappability File

Optional parameters:
Output File (Default based on BAM and BED file names)
Log2 output minimum threshold (Default -3)
Mappability threshold to consider (Default 0.25)

Example: python2 detect_deletion_BAM.py -b NGS.bam -c Interval.bed -m Mappability.tab -o NGS_deletion.tab -l -2 -M 0.25
"""

def calculateDeletion(PASS):
	SCORES = []
	for key in PASS:
		SCORES.append(float(PASS[key]))
	MEDIAN = np.nanmedian(SCORES)
	SCORES = []
	for key in PASS:
		if float(PASS[key]) == 0:
			SCORES.append((key, 'No Data Detected'))
                elif np.isnan(float(PASS[key])):
                        SCORES.append((key, 'Region does not meet mappability threshold'))
                else:
                        SCORES.append((key, np.log(float(PASS[key]) / MEDIAN) / np.log(2)))
	return SCORES

def iterateBAM(bam, bed, READLENGTH, MAP, MAPTHRESH):
	# open BAM file
	samfile = pysam.AlignmentFile(bam, "rb")
       	# open BED file
	file = open(bed, "r")

	# Counter of total tags mapping across all intervals
        totalTags = 0;
        totalSize = 0;

	PASS = {}
	FAIL = []
        totalTags = [0] * len(READLENGTH)

	# Iterate bed coord file, getting tag counts across interval
	for line in file:
		bedline = line.rstrip().split("\t")
		intervalID = bedline[3]
		intervalSize = int(bedline[2]) - int(bedline[1])
		if intervalSize < 1:
			FAIL.append(line.rstrip()+"|Interval coordinates invalid")
		try:
	                intervalCount = [0] * len(READLENGTH)
			mapAvg = [0] * len(READLENGTH)
			# Populate intervalCount array with tags in current BED region seperated by read length
			for read in samfile.fetch(bedline[0], int(bedline[1]), int(bedline[2])):
				if read.mapping_quality >= 5:
					index = closestLength(READLENGTH, read.query_length)
					totalTags[index] = totalTags[index] + 1
					intervalCount[index] = intervalCount[index] + 1
	                totalSize = totalSize + intervalSize

        	        # Calculate avg tags per bp across the entire region
	                intervalAvg = list(map(lambda x : float(x) / float(intervalSize), intervalCount))
			
	                # Normalize avg reads per interval by mappability
        	        for index in range(0, len(MAP[intervalID])):
                	        if float(MAP[intervalID][index]) >= MAPTHRESH:
                        	        mapAvg[index] = float(intervalAvg[index] * intervalCount[index]) / float(MAP[intervalID][index])
	                        else:   
        	                        mapAvg[index] = float('NaN')
			
			PASS[intervalID] = (mapAvg, intervalCount)

		except (IndexError, ValueError):
			FAIL.append(line.rstrip()+"|Not present in BAM File")

	SCORE = {}
	for key in PASS:
		# If read length was not detected in use across experiment, then consider it 'unmappable' and remove from consideration
		for index in range(0, len(MAP[key])):
			if float(totalTags[index]) == 0:
				MAP[key][index] = '0.0'
		# Calculate weighted average based on the variable read length
                if all(float(x) < MAPTHRESH for x in MAP[key]):
	                normalizedScore = float('NaN')
        	elif sum(PASS[key][1]) != 0:
                	normalizedScore = np.nansum(PASS[key][0]) / sum(PASS[key][1])
                else:
	                normalizedScore = 0
                SCORE[key] = normalizedScore


	if float(totalSize) <= 0:
		print "ERROR!!!\tTotal size of all intervals surveyed is less than 1"
		sys.exit(-1)
	
	# Close files
	file.close()
	samfile.close()
	
	return SCORE,FAIL

def closestLength(READLENGTH, read):
	index = -1
	dist = 99999
	for i in range(len(READLENGTH)):
		if(abs(read - int(READLENGTH[i])) < dist):
			dist = abs(read - int(READLENGTH[i]))
			index = i
	return index

def loadMap(MAP):
	# open Mappability file
        file = open(MAP, "r")
	header = 0;
	MAP = {}
        # Iterate BED coord file, getting tag counts across interval
        for line in file:
                mapline = line.rstrip().split("\t")
		if header == 0:
			header = mapline[1:]
		else:
			MAP[mapline[0]] = mapline[1:]
	file.close()
	return header,MAP

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

	BAM = BED = MAP = OUT = ""

        # Variable to set the mappability threshold so that we do not consider regions with mappability 
        # below this number 0-1, Default to 0.25 meaning at least 25% of the region must be uniquely mappable
        # by at least one actively used readlength
	MAPTHRESH = 0.25

	OUTPUTTHRESH = -3

	# get arguments
        optlist, alist = getopt.getopt(sys.argv[1:], 'hb:c:m:M:o:l:')
        for opt in optlist:
                if opt[0] == "-h": sys.exit(usage)
		elif opt[0] == "-b": BAM = opt[1]
                elif opt[0] == "-c": BED = opt[1]
		elif opt[0] == "-m": MAP = opt[1]
		elif opt[0] == "-M": MAPTHRESH = float(opt[1])
		elif opt[0] == "-o": OUT = opt[1]
		elif opt[0] == "-l": OUTPUTTHRESH = float(opt[1])
		else: sys.exit(usage)

        if BAM == "":
                print "No BAM file detected!!!"
                sys.exit(usage)
        elif BED == "":
                print "No BED Coordinate file detected!!!"
                sys.exit(usage)
        elif MAP == "":
                print "No Mappability file detected!!!"
                sys.exit(usage)
	if OUT == "":
		OUT = os.path.splitext(os.path.basename(BAM))[0] + "_" + os.path.splitext(os.path.basename(BED))[0] + ".tab"


        print "BAM file: ",BAM
        print "BED file: ",BED
        print "Mappability file: ",MAP
	print "Mappability threshold: ",MAPTHRESH
        print "Output file: ",OUT
	print "Log2 output threshold: ",OUTPUTTHRESH

        # Validate BAM file
        if(not validateBAM(BAM)):
                print "ERROR!!!\tNo BAM index detected.\n"
                sys.exit(-1)
	
	# Load mappability file
	READLENGTH, REGIONMAP = loadMap(MAP)
	print "Mappability file loaded"

	# Load BED file and calculate BAM coverage
	PASS, FAIL = iterateBAM(BAM, BED, READLENGTH, REGIONMAP, MAPTHRESH)
	print "Genomic coordinate coverage calculated"

	# Calculate log2 tag enrichment over median of mappability-normalized tag avg per region
	SCORE = calculateDeletion(PASS)	
	print "Depletion calculated"

	# Output final data
	output = open(OUT, "w")
        FINAL = sorted(SCORE, key=lambda x:x[1], reverse=True)
	for id,score in FINAL:
		if str(score) == 'No Data Detected':
			output.write(id + "\t" + str(score) + "\n")
        for id,score in reversed(FINAL):
                try:
			if float(score) < OUTPUTTHRESH:
                        	output.write(id + "\t" + str(score) + "\n")
			else:
				break
		except(ValueError):
			pass

	if FAIL:
		output.write("Invalid coordinates detected:\n")
		for line in FAIL:
			output.write(line + "\n")
		output.close()
