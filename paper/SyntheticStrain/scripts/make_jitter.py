from os import listdir
from os.path import isfile, join
import sys
import re
import random
import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Python 3.6+
# relies on dict insertion order

# Check Seaborn documentation: https://seaborn.pydata.org/generated/seaborn.swarmplot.html

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-i','--input-files', metavar='file_list', dest='file_list', required=True, help='Script takes text file with a list of file paths to the parsed results of each simulation experiment "<genome>_<strain>\t<depth>\t..." (i.e. "depth_simulations.txt")')
	parser.add_argument('-t','--title', metavar='figure_title', dest='title', required=True, help='Title to add to figure which is also used to name output image file')

	args = parser.parse_args()
	return(args)

def parse_data(data_file):
	'''Parse the data file with the simulation results (<genome>_<strain>_<depth>_scores.txt)'''
	data = {"value":[],"strain":[]}
	index_keys = []
	reader = open(data_file,'r')
	for line in reader:
		tokens = line.strip().split("\t")
		if(tokens[0].find("#")==0):
			# Initialize strain keys
			index_keys = tokens[1:]
			for strain in index_keys: data.update({strain:[]})
			continue
		for i in range(len(index_keys)):
			data["strain"].append(index_keys[i])
			# Set visualize-able values for "Inf" and "NaN"
			if(tokens[i+1]=="Inf"):
				data["value"].append(20.0)
			elif(tokens[i+1]=="NaN"):
				data["value"].append(-10.0)
			else:
				# Parse as float
				data["value"].append(float(tokens[i+1]))
	reader.close()

	return(data)

if __name__ == "__main__":
	'''Plot jitter/stripplot in R fashion using seaborn library'''
	args = getParams()
	# Initialize variables
	SIZE = 2
	JITTER = 1
	all_data = {"experiment":[],"value":[],"strain":[]}
	# Parse list of file paths
	i_reader = open(args.file_list,'r')
	for line in i_reader:
		data_file = line.strip().split("\t")[0]
		# Parse datafile for values
		data = parse_data(data_file)
		# Merge values into master dataframe
		all_data["experiment"].extend([data_file]*len(data["value"]))
		all_data["value"].extend(data["value"])
		all_data["strain"].extend(data["strain"])
	i_reader.close()

	# 1000 x 2 x 6 = 12000 points to plot
	# Plot data points using seaborn and label axes
	ax = sns.stripplot(x="experiment", y="value", hue="strain", data=all_data, size=SIZE, jitter=JITTER)
	ax.set_ylabel("log2 score")
	ax.set_xlabel("Simulation Experiment")
	ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
	# View/save plot
	#plt.show()
	out_png_fn = args.title.replace(" ","_")+".svg"
	plt.savefig(out_png_fn)
	#print(out_png_fn)
