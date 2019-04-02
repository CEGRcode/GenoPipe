import sys, getopt, os
from scipy.stats import poisson

usage = """
Usage:
This script takes in the paired-end table containing the stats from the EpitopeID pipeine and calculates
the Poisson p-value threshold given a minimum 2-fold enrichment over the expected tags based on the bin size.

Example: python2 calculate_EpitopeSignificance.py -t PE_table.out -p 0.05 -c 10000 -s 12000000 -o PE_sig.out
"""

# Main program which takes in input parameters
if __name__ == '__main__':
        if len(sys.argv) < 2 or not sys.argv[1].startswith("-"): sys.exit(usage)

	TABLE = PVAL = COUNT = SIZE = OUT = ""

	# get arguments
        optlist, alist = getopt.getopt(sys.argv[1:], 'ht:p:c:s:o:')
        for opt in optlist:
                if opt[0] == "-h": sys.exit(usage)
		elif opt[0] == "-t": TABLE = opt[1]
		elif opt[0] == "-p": PVAL = float(opt[1])
                elif opt[0] == "-c": COUNT = float(opt[1])
                elif opt[0] == "-s": SIZE = float(opt[1])
		elif opt[0] == "-o": OUT = opt[1]
		else: sys.exit(usage)

	if TABLE == "":
		print "No PE_table detected!!!"
		sys.exit(usage)
	elif PVAL == "":
		print "No Pvalue input!!!"
		sys.exit(usage)
        elif COUNT == "":
                print "No Single-end epitope counts input!!!"
                sys.exit(usage)
        elif SIZE == "":
                print "No Genome-size input!!!"
                sys.exit(usage)
	if OUT == "":
                OUT = os.path.splitext(os.path.basename(BAM))[0] + ".tab"

	# Minimum fold enrichment over background
	MINFOLD = 2;

	# Open output file for writing
        output = open(OUT, "w")
       	# open PE_table
	file = open(TABLE, "r")

	TABLE = []
	# Iterate table file, calculating pvalue as we go
	header = ""
	for line in file:
		if "detected" in line:
			output.write(line + "\n")
			file.close()
		elif header == "":
			header = "\nGeneID\tEpitopeID\tEpitopeLocation\tEpitopeCount\tpVal\n";
			output.write(header)
		if header != "":
			tableline = line.rstrip().split("\t")
			coord = tableline[0].split(":")
			WIDTH = coord[1].split("-")
			X = int(tableline[3])
			mu = (MINFOLD * (float(WIDTH[1]) - float(WIDTH[0])) * COUNT) / SIZE
			prob = 1 - poisson.cdf(X, mu) + poisson.pmf(X, mu)
			TABLE.append((line.rstrip(), prob))
	file.close()
	TABLE.sort(key=lambda tup: tup[1])
	for line,pval in TABLE:
		if float(pval) < PVAL:
			output.write(line + "\t" + str(pval) + "\n")
	output.close()
