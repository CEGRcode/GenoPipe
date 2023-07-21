#!/bin/python
from os import listdir
from os.path import isfile, join
import sys
import re
import random
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

# Python 3 needed for encoding feature for UTF-8
# (ENCODE uses some capital delta chars in summary descriptions of GeneticModifications)

# Check Seaborn documentation: https://seaborn.pydata.org/generated/seaborn.swarmplot.html

ENCODEtoStrainID = {
	"HeLa-S3":"HELA",
	"LNCAP":"LNCAPCLONEFGC",
	"MCF-7":"MCF7",
	"SK-N-SH":"SKNSH"
}

# K562       2640
# A549       1189
# MCF-7       556
# SK-N-SH     210
# HeLa-S3     196
# HCT116       96


def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Build histogram characterization plot from hg38 ENCODE StrainID results.')
	parser.add_argument('-m','--metadata', metavar='metadata_fn', required=True, help='the metadata file downloaded with ENCODE dataset')
	parser.add_argument('-i','--input-dir', metavar='input_dir', required=True, help='the directory where all the StrainID output files were saved (*strain.tab)')

	parser.add_argument('-o','--output', metavar='png_fn', required=True, help='the output figure image')

	parser.add_argument('-a','--assay', metavar='assay_name', default=None, help='the ENCODE assay name to filter datasets by (default:No Filter)')

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
def parse_file(var_file, expected):
	# Parse file
	scores = pd.read_table(var_file, sep='\t', header=0, names=['Strain','Scores'])
	# Add filename info
	scores['Filename'] = var_file
	# Add match information
	scores['Match'] = scores['Strain']==expected
	# Return scores
	return(scores)

if __name__ == "__main__":
	'''Plot swarm'''
	args = getParams()

	# Hardcode some presets
	SIZE=5
	strains2filter = ['P2URK562.vcf']
	strains2filter.extend(['HCT15.vcf', 'HCT8.vcf'])
	strains2filter.extend(['HEL9217.vcf', 'HEL.vcf'])
	strains2filter.extend(['HEP3B217.vcf'])
	strains2filter.extend(['MCF10A.vcf', 'MCF12A.vcf'])
	strains2filter.extend(['LN18.vcf', 'LN215.vcf', 'LN229.vcf', 'LN235.vcf', 'LN319.vcf', 'LN340.vcf', 'LN382.vcf', 'LN405.vcf', 'LN428.vcf', 'LN443.vcf', 'LN464.vcf', 'LNZTA3WT4.vcf', 'LNZ308.vcf'])

	strains2filter.extend(['SKN3.vcf', 'SKNAS.vcf', 'SKNBE2.vcf', 'SKNDZ.vcf', 'SKNEP1.vcf', 'SKNFI.vcf', 'SKNMC.vcf', 'SKNMM.vcf', 'SKNO1.vcf','SKN.vcf'])
	strains2filter.extend(['SKBR3.vcf', 'SKBR5.vcf', 'SKBR7.vcf', 'SKCO1.vcf', 'SKES1.vcf', 'SKGIIIA.vcf', 'SKGII.vcf', 'SKGI.vcf', 'SKGT2.vcf', 'SKGT4.vcf',
			'SKHEP1.vcf', 'SKLMS1.vcf', 'SKLU1.vcf', 'SKM1.vcf', 'SKMEL19.vcf', 'SKMEL1.vcf', 'SKMEL24.vcf', 'SKMEL28.vcf', 'SKMEL2.vcf', 'SKMEL30.vcf',
			'SKMEL31.vcf', 'SKMEL3.vcf', 'SKMEL5.vcf', 'SKMES1.vcf', 'SKMG1.vcf', 'SKMM2.vcf', 'SKOV3.vcf', 'SKPNDW.vcf', 'SKRC20.vcf', 'SKRC31.vcf', 'SKUT1.vcf'])


	# Parse metadata
	filedata = pd.read_csv(args.metadata, sep='\t', header=1)
	filedata['BIOSAMPLE_NAME'] = None
	df_list_scores = []

	# Loop through each sample
	for index, row in filedata.iterrows():
		# Map ENCODE-formatted strain to StrainID-formatted
		filedata['BIOSAMPLE_NAME'][index] = ENCODEtoStrainID.get(filedata['Biosample name'][index], filedata['Biosample name'][index])
		expected_vcfname = filedata['BIOSAMPLE_NAME'][index] + ".vcf"

		# Check file exists
		id_file = join(args.input_dir,"%s_strain.tab" % filedata['Accession'][index])
		if(not isfile(id_file)):
			continue

		# Parse ID file and add scores to final dataframe
		scores = parse_file(id_file, expected_vcfname)
		df_list_scores.append(scores)

		# for FCL in strains2filter:
		# 	print(scores[scores['Strain']==FCL])


	# Concatenate the strains together
	all_scores = pd.concat(df_list_scores)

	# Apply a hardcoded filter for parental strains

	# all_scores = all_scores.loc[all_scores['Strain'] in strains2filter]
	for FCL in strains2filter:
		print(FCL)
		# print(all_scores[all_scores['Strain']==FCL])
		all_scores = all_scores[all_scores['Strain']!=FCL]


	# Get counts for samples with all NaNs and format for table
	data_nans =  pd.DataFrame(all_scores[pd.isnull(all_scores['Scores'])])
	data_value =  pd.DataFrame(all_scores[~pd.isnull(all_scores['Scores'])])

	# print(data_nans['Strain'].value_counts())

	# Plot violin, swarms, and table
	fig, ax = plt.subplots()
	sns.histplot(ax=ax, x="Scores", binwidth=.1, data=data_value[~data_value['Match']], color='cyan')
	ax2 = ax.twinx()
	sns.histplot(ax=ax2, x="Scores", binwidth=.1, data=data_value[data_value['Match']], color='orange')
	plt.tight_layout()
	# palette = {
	# 	"A549":"tab:blue",
	# 	"HCT116":"tab:orange",
	# 	"HELA":"tab:green",
	# 	"HepG2":"tab:red",
	# 	"K562":"tab:purple",
	# 	"LnCap":"tab:olive",
	# 	"MCF7":"tab:cyan",
	# 	"SKnSH":"tab:pink"
	# }

	# Format figure
	ax.set_xlabel("StrainID -log2 score")
	ax.set_ylabel("number of scores for every sample x other cell lines (blue)")
	ax2.set_ylabel("number of scores for sample x matching cell line (orange)")
	# Save figure
	fig.set_size_inches(12,8)
	#plt.show()
	plt.savefig(args.output, dpi=500)
