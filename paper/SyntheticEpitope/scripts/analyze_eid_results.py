#!/bin/python
import os, sys, argparse

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Parse *-ID.tab files for epitope and localization information.')

	parser.add_argument('-i','--input-dir', metavar='input_dir', required=True, help='the directory where all the EpitopeID output files were saved (*-ID.tab)')
	parser.add_argument('-e','--epitope', metavar='expected_epitope', required=True, help='the expected epitope that should be found within the sequences')
	parser.add_argument('-t','--target', metavar='expected_target', required=True, help='the expected protein target the epitope should localize to')

	args = parser.parse_args()
	return(args)

'''
EpitopeID	EpitopeCount
Extended-Tap	14

GeneID	EpitopeID	EpitopeLocation	EpitopeCount	pVal
Epitope could not be detected genomically
'''

'''
EpitopeID	EpitopeCount
Extended-Tap	243
HA_v3	9

GeneID	EpitopeID	EpitopeLocation	EpitopeCount	pVal
HAT2|chr5:47168-48373	HA_v3	N-term	3	2.00151879412e-05
MHF1-Promoter|chr15:158923-159173	Extended-Tap	N/A	1	0.0103107852124
'''

# Main program which takes in input parameters
if __name__ == '__main__':
	'''Collect metadata and EpitopeID results to get detection stats on the YEP data'''

	args = getParams()
	nfiles = 0

	# Start copying metadata file with extra ID info appended
	for id_file in os.listdir(args.input_dir):
		# Initialize output vals
		found = 0
		localized = 0
		# Parse EpitopeID results files
		reader = open(args.input_dir + '/' + id_file,'r')
		line = "START"
		while( line!= "" ):
			line = reader.readline()
			tokens = line.strip().split('\t')
			# Parse epitope tag count info
			if (line.find("EpitopeID\tEpitopeCount")==0 ):
				while( len(tokens)==2 ):
					line = reader.readline()
					tokens = line.strip().split('\t')
					if(tokens[0]==args.epitope):
						found = tokens[1]
			# Parse localization details
			elif (line.find("GeneID\tEpitopeID\tEpitopeLocation")==0 ):
				while (len(tokens)==5):
					line = reader.readline()
					tokens = line.strip().split('\t')
					locus = tokens[0].strip().split("|")[0].split('-')[0].upper()
					if (locus==args.target.upper()):
						localized += 1
		reader.close()
		# Write results and update summary counts
		sys.stdout.write('%s\t%s\t%i\n' % (args.input_dir+id_file, found, localized ) )
		# Increment file count
		nfiles += 1
	reader.close()

	sys.stderr.write('Summary: %i files parsed.\n' % nfiles)
