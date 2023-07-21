import sys, os
import pandas as pd
import argparse

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Merge all tab files from StrainID results into a single table.')
	parser.add_argument('-m','--metadata', metavar='metadata_fn', required=True, help='the metadata file downloaded with ENCODE dataset')
	parser.add_argument('-i','--input-dir', metavar='input_dir', required=True, help='the directory where all the StrainID output files were saved (*strain.tab)')
	parser.add_argument('-o','--output', metavar='txt_fn', required=True, help='the output tab-delimited table path')

	args = parser.parse_args()
	return(args)

if __name__ == "__main__":
	'''Merge results'''

	args = getParams()

	merged_df = None
	# Loop through each StrainID results file
	reader = open(args.metadata, 'r')
	for line in reader:
		if(line.find('/files/')!=0):
			continue
		tokens = line.split('\t')
		# Load results and name VCF filename column "CellLine"
		temp_df = pd.read_csv(args.input_dir + "/" + tokens[1] + "_strain.tab", sep='\t', header=0)
		temp_df.columns.values[0] = 'CellLine'
		temp_df.loc[len(temp_df.index)] = ['0Correct_Strain', tokens[12]]
		# Merge into merged_df
		if (merged_df is None):
			merged_df = temp_df
			continue
		merged_df = pd.merge(merged_df, temp_df, on='CellLine')
	reader.close()
	# sort values before saving
	merged_df = merged_df.sort_values(by=['CellLine'])
	# write merged_df to output file
	merged_df.to_csv(args.output, sep='\t')
