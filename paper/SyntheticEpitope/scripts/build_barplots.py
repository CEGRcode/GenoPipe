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
#       Depth Epitope  Target
# 0      100K    R500    CTCF
# 1      100K    R100    CTCF
# 2      100K     R50  POLR2H
# 3      100K     R20    CTCF
# 4      100K     R20  POLR2H
# ...     ...     ...     ...
# 30360   50M     R50     YY1
# 30361   50M     R50     YY1
# 30362   50M     R50     YY1
# 30363   50M     R50     YY1
# 30364   50M     R50     YY1

	localized = pd.DataFrame(np.where(filedata['LocalizationCount']==0, False, True), columns=['Localized'])
#        Localized
# 0           True
# 1           True
# 2           True
# 3           True
# 4           True
# ...          ...
# 30360       True
# 30361      False
# 30362       True
# 30363       True
# 30364       True

	bardata = pd.concat([filedata,experiment_info,localized],axis=1)
	print(bardata)
	filedata = experiment_info = localized = None

	# Filter down to localized data counts
	localized_bardata = bardata[bardata['Localized']]
	print(localized_bardata)
	# Reorganize
	localized_counts = localized_bardata.groupby(['Target','Epitope','Depth']).size().reset_index().rename(columns={0:'count'})
	print(localized_bardata)
	# print(localized_counts)

	# Hardcode depth order list for yeast or human simulations
	depth_order = ["1M","10M","20M","50M"]
	if (args.yeast):
		depth_order = ["10K","100K","1M","10M"]
	target_order = ["CTCF","POLR2H","YY1"]
	if (args.yeast):
		target_order = ["Rap1","Reb1","Sua7"]

	# Configure and plot data into grouped bars
	pal = sns.color_palette("viridis", len(pd.unique(localized_counts['Target'])))
	if (args.yeast):
		pal = sns.color_palette("YlOrBr", len(pd.unique(localized_counts['Target'])))

	# Super plot
	fig, axes = plt.subplots(2,2, sharex=True, sharey=True)
	if (args.yeast):
		fig.suptitle('Yeast')
	else:
		fig.suptitle('Human')

	# Plot subplots by Synthetic epitope length
	for i,ename in enumerate(["R500", "R100", "R50", "R20"]):
		axes[i//2,i%2].set_title(ename)
		subplotdata = localized_counts[localized_counts['Epitope'] == ename]
		print(subplotdata)
		# print(subplotcounts)
		sns.barplot(ax=axes[i//2,i%2], data=subplotdata, x='Depth', y='count', hue='Target', hue_order=target_order, order=depth_order, palette=pal, ci=None, bottom=0)
	# fig = plot.get_figure()
	fig.savefig(args.output)
