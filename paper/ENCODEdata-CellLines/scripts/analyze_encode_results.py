from os import listdir
from os.path import isfile, join
import sys
import argparse
import numpy as np
import pandas as pd

# Python 3 needed for encoding feature for UTF-8
# (ENCODE uses some capital delta chars in summary descriptions of GeneticModifications)

CL = ["LnCap", "MCF7", "SKnSH", "HepG2", "K562", "A549", "HCT116", "HELA"]
ENCODEtoStrainID = {
	"HeLa-S3":"HELA",
	"LNCAP":"LnCap",
	"MCF-7":"MCF7",
	"SK-N-SH":"SKnSH"
}

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Parse metadata file and StrainID output to check per sample detection by StrainID scores.')
	parser.add_argument('-m','--metadata', metavar='metadata_fn', required=True, help='the metadata file downloaded with ENCODE dataset that includes info like PE/SE, cell line, assay type, and read lengths/SE-PE')
	parser.add_argument('-i','--input-dir', metavar='input_dir', required=True, help='the directory where all the EpitopeID output files were saved (*strain.tab)')
	parser.add_argument('-o','--output', metavar='output_fn', required=True, help='the output filepath for final TSV with parsed StrainID scores')

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
	scores = []
	reader = open(var_file,'r')
	for line in reader:
		tokens = line.split("\t")
		if(tokens[0]==""):
			continue
		score = float(tokens[1].strip())
		if(tokens[1].strip().lower()=="inf"):
			score = np.Inf
		elif(tokens[1].strip().lower()=="nan"):
			score = np.NaN
		# update dict
		scores.append((score, tokens[0].split(".")[0]))
	reader.close()
	return(scores)

#ENCFF364CPX	/files/ENCFF364CPX/@@download/ENCFF364CPX.bam	79651f67b1c4d564395c18be9cdff62f	HeLa-S3	ChIP-seq	single-ended	36	/files/ENCFF807MUK/|/files/ENCFF000BAO/	unfiltered alignments	/experiments/ENCSR000AOB/	/biosample-types/cell_line_EFO_0002791/	1884977463	released	2020-02-18T20:47:36.519163+00:00
#ENCFF325UJS	/files/ENCFF325UJS/@@download/ENCFF325UJS.bam	042b20b3e149df6c1f4e5c95f83653ee	HepG2	ChIP-seq	single-ended	36	/files/ENCFF807MUK/|/files/ENCFF000BGR/	alignments	/experiments/ENCSR000AOM/	/biosample-types/cell_line_EFO_0001187/	978259147	released	2020-02-18T09:15:21.891603+00:00
#ENCFF821WQW	/files/ENCFF821WQW/@@download/ENCFF821WQW.bam	e9c6eeedee7dc41d6e19ca6f7a6777f3	HepG2	ChIP-seq	paired-ended	100	/files/ENCFF195VBJ/|/files/ENCFF807MUK/|/files/ENCFF594PZU/	alignments	/experiments/ENCSR730TBC/	/biosample-types/cell_line_EFO_0001187/	4831825733	released	2017-06-06T18:20:14.545786+00:00
if __name__ == "__main__":
	'''Collect metadata and StrainID results to get detection stats on the cell line ENCODE data'''
	args = getParams()

	# Parse metadata
	data = pd.read_csv(args.metadata, sep='\t', names=['File_Accession','Download_URL','MD5sum', 'ENCODE_strain', 'Assay', 'library_type', 'read_length', 'derived_from', 'bam_type', 'Experiment_Accession', 'Biosample_Accession', 'File_Size', 'Audit_status', 'date'])

	# Initialize summary StrainID results columnss
	data['StrainID_strain'] = None
	data['StrainID_success'] = None
	data['Comment'] = None

	# Initialize each strain's score to Nones
	for c in CL:
		data[c + "_score"] = None

	# Loop through each sample
	for index, row in data.iterrows():
		# Map ENCODE-formatted strain to StrainID-formatted
		data['ENCODE_strain'][index] = ENCODEtoStrainID.get(data['ENCODE_strain'][index], data['ENCODE_strain'][index])
		# Check file exists
		id_file = join(args.input_dir,"%s_strain.tab" % data['File_Accession'][index])
		# print(data['File_Accession'][index])
		if(not isfile(id_file)):
			data['StrainID_success'][index] = "Missing Results"
			data['Comment'][index] = "Missing Results"
			continue
		# Parse id file for cell line score info and sort strains by score
		strain_info = parse_file(id_file)
		for s in strain_info:
			data[s[1] + "_score"][index] = str(s[0])
		# Sort scores to get the best one
		strain_sortbyscore = sorted([s for s in strain_info if not np.isnan(s[0]) ], reverse=True)
		# Assign strain with best score to called_strain
		if(len(strain_sortbyscore)>0):
			data['StrainID_strain'][index] = strain_sortbyscore[0][1]
			data['StrainID_success'][index] = data['StrainID_strain'][index] == data['ENCODE_strain'][index]

	# Write final data frame
	data.to_csv(args.output, sep="\t")
