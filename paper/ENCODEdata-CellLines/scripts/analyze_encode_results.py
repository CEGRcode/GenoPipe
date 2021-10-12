from os import listdir
from os.path import isfile, join
import sys
import argparse

# Python 3 needed for encoding feature for UTF-8
# (ENCODE uses some capital delta chars in summary descriptions of GeneticModifications)

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Parse metadata file and GenoPipe output to check detection rates of the GenoPipe tool.')
	parser.add_argument('-m','--metadata', metavar='metadata_fn', required=True, help='the metadata file downloaded with ENCODE dataset that includes info like PE/SE, cell line, assay type, and read lengths/SE-PE')
	parser.add_argument('-i','--input-dir', metavar='input_dir', required=True, help='the directory where all the EpitopeID output files were saved (*strain.tab)')

	args = parser.parse_args()
	return(args)

#	ENCFF000DZC.bam
#LnCap.vcf	-5.082117812158647
#MCF7.vcf	-6.1143012059424935
#SKnSH.vcf	-5.7641601741217645
#HepG2.vcf	-5.595833186705702
#K562.vcf	1.8812984639660986
#A549.vcf	-6.059318584695944
#HCT116.vcf	-4.847450343904915
#HELA.vcf	-4.906670711358038
def parse_file(var_file):
	dict = {}
	reader = open(var_file,'r')
	for line in reader:
		tokens = line.split("\t")
		if(tokens[0]==""):
			continue
		score = float(tokens[1].strip())
		if(tokens[1].strip().lower()=="inf"):
			score = 500000
		elif(tokens[1].strip().lower()=="nan"):
			score = -500000

		# update dict
		dict[tokens[0].split(".")[0]] = score
	reader.close()
	return(dict)

#ENCFF364CPX	/files/ENCFF364CPX/@@download/ENCFF364CPX.bam	79651f67b1c4d564395c18be9cdff62f	HeLa-S3	ChIP-seq	single-ended	36	/files/ENCFF807MUK/|/files/ENCFF000BAO/	unfiltered alignments	/experiments/ENCSR000AOB/	/biosample-types/cell_line_EFO_0002791/	1884977463	released	2020-02-18T20:47:36.519163+00:00
#ENCFF325UJS	/files/ENCFF325UJS/@@download/ENCFF325UJS.bam	042b20b3e149df6c1f4e5c95f83653ee	HepG2	ChIP-seq	single-ended	36	/files/ENCFF807MUK/|/files/ENCFF000BGR/	alignments	/experiments/ENCSR000AOM/	/biosample-types/cell_line_EFO_0001187/	978259147	released	2020-02-18T09:15:21.891603+00:00
#ENCFF821WQW	/files/ENCFF821WQW/@@download/ENCFF821WQW.bam	e9c6eeedee7dc41d6e19ca6f7a6777f3	HepG2	ChIP-seq	paired-ended	100	/files/ENCFF195VBJ/|/files/ENCFF807MUK/|/files/ENCFF594PZU/	alignments	/experiments/ENCSR730TBC/	/biosample-types/cell_line_EFO_0001187/	4831825733	released	2017-06-06T18:20:14.545786+00:00
if __name__ == "__main__":
	'''Collect metadata and StrainID results to get detection stats on the cell line ENCODE data'''
	args = getParams()

	# Parse metadata
	reader = open(args.metadata, 'r', encoding='utf-8')
	for mline in reader:
		# Pull relevant info from metadata tokens
		mtokens = mline.strip().split('\t')
		encff = mtokens[0]

		# Check file exists
		id_file = join(args.input_dir,"%s_strain.tab" % encff)
		if(not isfile(id_file)):
			sys.stderr.write("%s: no results generated.\n" % (id_file))
			continue

		# Initialize id file variables to save
		called_strain = ""
		# Parse id file for cell line score info and sort strains by score
		strain_info = parse_file(id_file)
		strain_sortbyscore = sorted( strain_info.keys(), key=lambda x: (strain_info[x]), reverse=True)
		# Assign strain with best score to called_strain
		if(len(strain_info.keys())>1):
			called_strain = strain_sortbyscore[0]

		# Write called strain with metadata
		sys.stdout.write( "%s\t%s\t%s\t%s\n" % (encff, called_strain, strain_info.get(called_strain,"NaN"), "\t".join(mtokens[3:])) )
	reader.close()
