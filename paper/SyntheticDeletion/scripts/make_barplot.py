from os import listdir
from os.path import isfile, join
import sys
import re
import argparse
import matplotlib.pyplot as plt
import numpy as np

# Python 3.6+
# relies on dict insertion order

# Check Matplotlib colors when building your config files: https://matplotlib.org/stable/gallery/color/named_colors.html

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='')
		
	parser.add_argument('-t','--title', metavar='figure_title', dest='title', required=True, help='')
	parser.add_argument('-c','--config', metavar='config_fn', dest='config_fn', required=True, help='the config file for a grid-organized subplot')
	parser.add_argument('-d','--header', dest='header', default=False, required=False, help='skip first line as column header')
	
	args = parser.parse_args()
	return(args)


def get_line_count(simulation_summary_file):
	line_count = 0
	reader = open(simulation_summary_file,'r')
	for line in reader:
		line_count += 1
	reader.close()
	return(line_count)

def parse_configs(configs_fn):
	#DEPTH_ORDER	 500K	600K	700K	800K	900K	1M
	#STRAIN_COLOR	Rap1-del		orange
	#STRAIN_COLOR	Reb1-del		yellow
	#filenameA	  Rap1-del		500K
	#filenameB	Reb1-del		1M
	# Dictionary Structure:
	# { "strainA":["color", [ (depth1,val1), (depth2,val2), ... ] ], "strainB": ...}
	subplot_configs = {}
	depth_order = []
	reader = open(configs_fn,'r')
	for line in reader:
		tokens = line.strip().split('\t')
		if(line.find("DEPTH_ORDER")==0):
			depth_order = tokens[1:]
			continue
		
		#subplot_configs.setdefault(tokens[1],["black",{}])
		subplot_configs.setdefault(tokens[1],["black",[]])
		if(line.find("STRAIN_COLOR")==0):
			subplot_configs[tokens[1]][0] = tokens[2]
			continue
		# GET VALUE (change to parse file later)
		#tokens[0] = get_line_count(tokens[0])
		tokens[0] = get_line_count(tokens[0])
		#subplot_configs[tokens[1]][1].update({tokens[2]:tokens[0]})
		subplot_configs[tokens[1]][1].append((tokens[2],tokens[0]))
	reader.close()
	

	#if(depth_order!=[]):
	#	depth_order = [  ]

	# Sort data according to DEPTH_ORDER
	for strain in subplot_configs.keys():
		data_list = subplot_configs[strain][1]
		subplot_configs[strain][1] = sorted(subplot_configs[strain][1], key = lambda x:depth_order.index(x[0]))
		#data = list(strain[1])
		#data_d = [ d[0] for d in data ]
		#data = [  for d in depth_order ]

	return(subplot_configs)

if __name__ == "__main__":
	'''Plot scatter'''
	args = getParams()
	BAR_WIDTH  = 0.25
	
	CONFIGS = parse_configs(args.config_fn)
	print(CONFIGS)
	
	fig = plt.figure()
	strain_list = list(CONFIGS.keys())
	for s in range(len(strain_list)):
		strain = strain_list[s]
		data = CONFIGS[strain][1]
		X = np.arange(len(data))
		plt.bar(X + s*BAR_WIDTH, [d[1] for d in data], color = CONFIGS[strain][0], width = BAR_WIDTH)
	plt.xticks(X, [ x[0] for x in CONFIGS[strain_list[0]][1]])
	plt.ylim(top=1000,bottom=0)
	plt.ylabel("Number of successful simulations\n(out of 1000)")
	plt.xlabel("Dataset size\n(# of paired-end reads sampled)")
	plt.title(args.title)
	#plt.show()
	out_png_fn = args.title.replace(" ","_")+".png"
	plt.savefig(out_png_fn)
	print(out_png_fn)
