#!/bin/python
import os, sys, argparse
import pandas as pd
import numpy as np
import seaborn as sns

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Build barplots from simulation results.')

	parser.add_argument('-i','--input', metavar='input_fn', required=True, help='the concatenated input from analyzed_eid_results.py')
	parser.add_argument('-o','--output', metavar='png_fn', required=True, help='the output figure image')

	parser.add_argument('-y','--yeast', action='store_true', help='flag to use yeast depths for x-axis')

	args = parser.parse_args()
	return(args)

'''
results/sacCer3/1M/R500/Reb1/ID//Simulation_953_R1-ID.tab	   77	  3
results/sacCer3/1M/R500/Reb1/ID//Simulation_668_R1-ID.tab	   89	  3
results/sacCer3/1M/R500/Reb1/ID//Simulation_923_R1-ID.tab	   73	  3
results/sacCer3/1M/R500/Reb1/ID//Simulation_641_R1-ID.tab	   73	  3
results/sacCer3/1M/R500/Reb1/ID//Simulation_581_R1-ID.tab	   89	  1
results/sacCer3/10K/R500/Reb1/ID/Simulation_5_R1-ID.tab 1	   1
results/sacCer3/10K/R500/Reb1/ID/Simulation_487_R1-ID.tab	   3	   0
results/sacCer3/10K/R500/Reb1/ID/Simulation_109_R1-ID.tab	   2	   0
results/sacCer3/10K/R500/Reb1/ID/Simulation_965_R1-ID.tab	   2	   0
results/sacCer3/10K/R500/Reb1/ID/Simulation_46_R1-ID.tab	0	   0
results/sacCer3/10K/R500/Reb1/ID/Simulation_49_R1-ID.tab	1	   1
results/sacCer3/10K/R500/Reb1/ID/Simulation_259_R1-ID.tab	   0	   0
results/sacCer3/10K/R500/Reb1/ID/Simulation_582_R1-ID.tab	   0	   0
results/sacCer3/10K/R500/Reb1/ID/Simulation_609_R1-ID.tab	   3	   1
'''

# Main program which takes in input parameters
if __name__ == '__main__':
	'''Collect metadata and EpitopeID results to get detection stats on the YEP data'''

	args = getParams()
	# bardata = pd.DataFrame({
	# 	'Filename':,
	# 	'Target':[],
	# 	'Epitope':[],
	# 	'Depth':[],
	# 	'EpitopeCount':[],
	# 	'LocalizationCount':[],
	# 	'Localized':[],
	# })

	# Populate dataframe with tab file data
	filedata = pd.read_table(args.input, sep='\t', names=['Filename','EpitopeCount','LocalizationCount'])
	experiment_info = filedata['Filename'].str.split('/',n=5,expand=True).loc[:,2:4]
	experiment_info = experiment_info.rename(columns={2:'Depth', 3:'Epitope', 4:'Target'})
	localized = pd.DataFrame(np.where(filedata['LocalizationCount']==0, False, True), columns=['Localized'])
	bardata = pd.concat([filedata,experiment_info,localized],axis=1)
	filedata = experiment_info = localized = None

	# Filter down to localized data counts
	localized_bardata = bardata[bardata['Localized']]
	localized_counts = localized_bardata.groupby(['Target','Epitope','Depth']).size().reset_index().rename(columns={0:'count'})
	print(localized_counts)

	# Hardcode depth order list for yeast or human simulations
	depth_order = ["100K","1M","10M","20M","50M"]
	if (args.yeast):
		depth_order = ["10K","100K","1M","10M"]

	# Configure and plot data into grouped bars
	pal = sns.color_palette("viridis", len(pd.unique(localized_counts['Target'])))
	# pal = sns.color_palette("YlOrBr", len(pd.unique(localized_counts['Target'])))
	# pal = sns.color_palette("Paired")
	plot = sns.barplot(data=localized_counts, x='Depth', y='count', hue='Target', order=depth_order, palette=pal, ci=None)
	fig = plot.get_figure()
	fig.savefig(args.output)
