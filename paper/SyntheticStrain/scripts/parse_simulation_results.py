from os import listdir
from os.path import isfile, join
import sys
import argparse

# Python 3 needed for encoding feature for UTF-8
# (ENCODE uses some capital delta chars in summary descriptions of GeneticModifications)

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Parse metadata file and GenoPipe output to check detection rates of the GenoPipe tool.')
	parser.add_argument('-i','--input-dir', metavar='input_dir', required=True, help='the directory where all the StrainID output files were saved (*strain.tab)')
	parser.add_argument('-v','--vcf-dir', metavar='vcf_dir', required=True, help='the directory where all the StrainID VCF db files are housed (for header formatting purposes)')
	args = parser.parse_args()
	return(args)

# 	Simulation_11.bam
# Y55.gatk.vcf	5.893664745278011
# BY4741.gatk.vcf	4.3921009156813495
# SEY6210.gatk.vcf	6.799588141350312
# Sigma1278b-10560-6B.gatk.vcf	6.236837493072714
# CEN.PK2-1Ca.gatk.vcf	8.360122378850173
# D273-10B.gatk.vcf	6.24252455963084
# RM11-1A.gatk.vcf	5.897733174105499
# BY4742.gatk.vcf	4.57267316132317
# FL100.gatk.vcf	6.480670542594632
# X2180-1A.gatk.vcf	4.525367446544814
# JK9-3d.gatk.vcf	6.523429642303322
# W303.gatk.vcf	6.593092702562104
def parse_file(var_file):
	dict = {}
	reader = open(var_file,'r')
	for line in reader:
		tokens = line.strip().split("\t")
		if(len(tokens)==1):
			continue
		dict[tokens[0].split(".")[0]] = tokens[1]
	reader.close()
	return(dict)


if __name__ == "__main__":
	'''Collect metadata and StrainID results to get detection stats on the cell line ENCODE data'''
	args = getParams()

	# Parse strains to track
	strain_keys = []
	for filename in listdir(args.vcf_dir):
		filename_tokens = filename.split(".")
		if(filename_tokens[-1]=="vcf"):
			strain_keys.append(filename_tokens[0])
	strain_keys.sort()

	# Write header
	sys.stdout.write("#\t%s\n" % "\t".join(strain_keys))

	# Parse metadata
	for sindex in range(1,1001):
		# Check file exists
		id_file = join(args.input_dir,"Simulation_%i_strain.tab" % sindex)
		if(not isfile(id_file)):
			sys.stderr.write("%s: no results generated.\n" % (id_file))
			continue

		# Parse id file for cell line score info
		strain_info = parse_file(id_file)

		# Write called strain with metadata
		sys.stdout.write( "Simulation_%i\t%s\n" % (sindex, "\t".join([strain_info.get(s,"-") for s in strain_keys]) ))
