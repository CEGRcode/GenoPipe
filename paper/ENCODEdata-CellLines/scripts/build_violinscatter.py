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

# Python 3.6+
# relies on dict insertion order

# Check Seaborn documentation: https://seaborn.pydata.org/generated/seaborn.swarmplot.html

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Build violinplots from ENCODE StrainID results.')

	parser.add_argument('-i','--input', metavar='input_fn', required=True, help='the output tab file from analyzed_eid_results.py')
	parser.add_argument('-o','--output', metavar='png_fn', required=True, help='the output figure image')

	parser.add_argument('-a','--assay', metavar='assay_name', default=None, help='the ENCODE assay name to filter datasets by (default:No Filter)')

	args = parser.parse_args()
	return(args)

if __name__ == "__main__":
	'''Plot swarm'''
	args = getParams()

	# Hardcode some presets
	SIZE=5
	palette = {
		"A549":"tab:blue",
		"HCT116":"tab:orange",
		"HELA":"tab:green",
		"HepG2":"tab:red",
		"K562":"tab:purple",
		"LnCap":"tab:olive",
		"MCF7":"tab:cyan",
		"SKnSH":"tab:pink"
	}

	# Parse data table results and get max StrainID scores
	filedata = pd.read_table(args.input, sep='\t')
	justscores = filedata.loc[:, 'LnCap_score':'HELA_score']
	filedata['StrainID_bestscore'] = justscores.max(axis=1)
	filedata = filedata.sort_values(by='ENCODE_strain')

	# Filter by assay name if specified
	if (args.assay!=None):
		filedata = filedata[filedata['Assay']==args.assay]

	# Separate success/fail sets for violin vs swarms
	data_success = filedata[filedata['StrainID_success']=='True']
	data_fails = filedata[filedata['StrainID_success']=='False']

	# Get counts for samples with all NaNs and format for table
	data_nans =  pd.DataFrame(filedata[filedata['StrainID_success'].isnull()]['ENCODE_strain'].value_counts(), index=["A549", "HCT116", "HELA", "HepG2", "K562", "LnCap", "MCF7", "SKnSH"]).T.fillna(value=0)

	# Plot violin, swarms, and table
	fig, ax = plt.subplots()
	sns.violinplot(ax=ax, x="ENCODE_strain", y="StrainID_bestscore", hue="StrainID_strain", data=data_success, palette=palette)
	sns.stripplot(ax=ax, x="ENCODE_strain", y="StrainID_bestscore", hue="StrainID_strain", data=data_fails, palette=palette)
	plt.table(cellText=data_nans.values, rowLabels=data_nans.index, colLabels=data_nans.columns)

	# Format figure
	ax.set_ylabel("StrainID -log2 score")
	ax.set_ylim(-8,10)

	# Save figure
	fig.set_size_inches(12,8)
	#plt.show()
	plt.savefig(args.output, dpi=500)
