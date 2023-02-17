#!/bin/python
import os, sys, argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Build barplots from simulation results.')

	parser.add_argument('-i','--input', metavar='input_fn', required=True, help='the concatenated input from analyzed_eid_results.py')
	parser.add_argument('-o','--output', metavar='png_fn', required=True, help='the output figure image')

	parser.add_argument('-y','--yeast', action='store_true', help='flag to use yeast depths for x-axis')

	args = parser.parse_args()
	return(args)

'''
results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_1000_R1-ID.tab	6	0	results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_1000_R1-ID.tab	6	0
results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_100_R1-ID.tab	20	3	results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_100_R1-ID.tab	20	0
results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_101_R1-ID.tab	10	2	results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_101_R1-ID.tab	10	0
results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_102_R1-ID.tab	4	1	results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_102_R1-ID.tab	4	0
results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_103_R1-ID.tab	9	2	results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_103_R1-ID.tab	9	0
results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_104_R1-ID.tab	7	3	results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_104_R1-ID.tab	7	0
results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_105_R1-ID.tab	6	2	results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_105_R1-ID.tab	6	1
results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_106_R1-ID.tab	9	2	results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_106_R1-ID.tab	9	1
results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_107_R1-ID.tab	10	2	results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_107_R1-ID.tab	10	0
results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_108_R1-ID.tab	12	3	results/mix_yeast/100K/90_Reb1_Rap1/ID/Simulation_108_R1-ID.tab	12	0
'''

# Main program which takes in input parameters
if __name__ == '__main__':
	'''Collect metadata and EpitopeID results to get detection stats on the YEP data'''

	args = getParams()

	# Populate dataframe with tab file data
	filedata = pd.read_table(args.input, sep='\t', names=['A_Filename','A_EpitopeCount','A_LocalizationCount','B_Filename','B_EpitopeCount','B_LocalizationCount'])
	# filedata = pd.read_table('results/MixSummaryReport_sacCer3_1M.txt', sep='\t', names=['A_Filename','A_EpitopeCount','A_LocalizationCount','B_Filename','B_EpitopeCount','B_LocalizationCount'])
	# Confirm parallel matching Filenames
	if (not filedata['A_Filename'].equals(filedata['B_Filename'])):
		print("Filename columns (A v B) do not match!!!")
		quit()

	# Drop unused fields
	filedata = filedata.drop(['B_Filename', 'A_EpitopeCount','B_EpitopeCount'], axis=1)

	# Parse out ratio info
	filedata['MixRatio'] = pd.to_numeric(filedata['A_Filename'].str.split('/',n=5,expand=True).loc[:,3].str.split('_',n=3,expand=True).loc[:,0])

	# Parse out Localization Counts
	filedata['A_Localized'] = np.where(filedata['A_LocalizationCount']==0, False, True)
	filedata['B_Localized'] = np.where(filedata['B_LocalizationCount']==0, False, True)
	filedata['Both_Localized'] = pd.concat([filedata['A_Localized'], filedata['B_Localized']], axis=1).all(axis=1)
	print(filedata)

	# Reorganize to long form
	tempA = filedata.groupby(['MixRatio','A_Localized']).size().reset_index().rename(columns={0:'count'})
	tempA['Type'] = 'A_Localized'
	tempA = tempA.rename({'A_Localized':'Localized'}, axis=1)
	tempB = filedata.groupby(['MixRatio','B_Localized']).size().reset_index().rename(columns={0:'count'})
	tempB['Type'] = 'B_Localized'
	tempB = tempB.rename({'B_Localized':'Localized'}, axis=1)
	tempBoth = filedata.groupby(['MixRatio','Both_Localized']).size().reset_index().rename(columns={0:'count'})
	tempBoth['Type'] = 'Both_Localized'
	tempBoth = tempBoth.rename({'Both_Localized':'Localized'}, axis=1)
	localized_linedata = pd.concat([tempA,tempB,tempBoth], ignore_index=True)
	print(localized_linedata)

	# Filter for localized counts
	localized_linedata = localized_linedata.drop(localized_linedata[localized_linedata['Localized']==False].index)
	filedata = tempA = tempB = tempBoth = None

	#   MixRatio  Localized  count            Type
	# 1       50       True    830     A_Localized
	# 3       60       True    866     A_Localized
	# 5       70       True    905     A_Localized
	# 7       80       True    935     A_Localized
	# 9       90       True    953     A_Localized
	# 1       50       True    835     B_Localized
	# 3       60       True    734     B_Localized
	# 5       70       True    654     B_Localized
	# 7       80       True    501     B_Localized
	# 9       90       True    294     B_Localized
	# 1       50       True    697  Both_Localized
	# 3       60       True    639  Both_Localized
	# 5       70       True    589  Both_Localized
	# 7       80       True    464  Both_Localized
	# 9       90       True    278  Both_Localized

	# Configure and plot data into grouped bars
	pal = sns.color_palette("viridis", 3)
	if (args.yeast):
		pal = sns.color_palette("YlOrBr", 3)
	pal[2] = "black"

	# Make Line plot
	snsplot = sns.lineplot(data=localized_linedata, x='MixRatio', y='count', hue='Type', hue_order=['A_Localized','B_Localized','Both_Localized'], palette=pal)#, ci=None, bottom=0)

	# Title the plot
	fig = snsplot.get_figure()
	if (args.yeast):
		fig.suptitle('Mix-Yeast')
	else:
		fig.suptitle('Mix-Human')

	# Set axis limits
	plt.xlim(0, 100)
	plt.ylim(0, 1000)
	plt.xticks(range(0,101,10))

	# Save figure by output filename
	# fig.savefig('results/blah.svg')
	fig.savefig(args.output)
