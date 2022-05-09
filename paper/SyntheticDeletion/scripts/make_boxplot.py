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


def parse_times(simulation_runtime_file):
#logs/depth.did.Rap1.500K.log.err-1:real    0m24.981s
#logs/depth.did.Rap1.500K.log.err-10:real    0m26.424s
#logs/depth.did.Rap1.500K.log.err-11:real    0m25.349s
#logs/depth.did.Rap1.500K.log.err-12:real    0m24.692s
	runtimes = []
	reader = open(simulation_runtime_file,'r')
	for line in reader:
		tokens = line.strip().split('\t')
		minutes = float(tokens[1].split('m')[0])
		seconds = float(tokens[1].split('m')[1].split('s')[0])
		runtimes.append(seconds + (60.0*minutes))
	reader.close()
	return(runtimes)

def parse_configs(configs_fn):
#DEPTH_ORDER     500K    600K    700K    800K    900K    1M
#STRAIN_ORDER   Rap1-del    Reb1-del
#STRAIN_COLOR    Rap1-del        orange
#STRAIN_COLOR    Reb1-del        yellow
#25      Rap1-del        500K
#1000    Reb1-del        1M
# Dictionary Structure:
# { "strainA":["color", [ (depth1,val1), (depth2,val2), ... ] ], "strainB": ...}
	subplot_configs = {}
	depth_order = []
	reader = open(configs_fn,'r')
	for line in reader:
		tokens = line.strip().split('\t')
		if(line.find("DEPTH_ORDER")==0):
			depth_order = tokens[1:]
			global n_depth
			n_depth = len(depth_order)
			continue
		#subplot_configs.setdefault(tokens[1],["black",{}])
		subplot_configs.setdefault(tokens[1],["black",[]])
		if(line.find("STRAIN_COLOR")==0):
			subplot_configs[tokens[1]][0] = tokens[2]
			continue
		# GET VALUE (change to parse file later)
		#tokens[0] = get_line_count(tokens[0])
		tokens[0] = parse_times(tokens[0])
		#subplot_configs[tokens[1]][1].update({tokens[2]:tokens[0]})
		subplot_configs[tokens[1]][1].append((tokens[2],tokens[0]))
	reader.close()
	
	#if(depth_order==[]):
		#depth_order = [ ]
	
	# Sort data according to DEPTH_ORDER
	n_strain = len(subplot_configs.keys())
	for strain in subplot_configs.keys():
		data_list = subplot_configs[strain][1]
		subplot_configs[strain][1] = sorted(subplot_configs[strain][1], key = lambda x:depth_order.index(x[0]))
	
	return(subplot_configs)

if __name__ == "__main__":
	'''Plot scatter'''
	args = getParams()
	
	CONFIGS = parse_configs(args.config_fn)
	#print(CONFIGS)
	BOX_WIDTH = 0.6
	
	global n_depth
	xticks = []
	xlabels = []
	strain_list = list(CONFIGS.keys())
	n_strains = len(strain_list)
	base_pos = range(0,n_depth*(n_strains+1),n_strains+1)
	
	fig = plt.figure()
	
	for s in range(n_strains):
		SCONFIG = CONFIGS[strain_list[s]]
		if(xlabels==[]):
			xlabels = [ SCONFIG[1][d][0] for d in range(len(SCONFIG[1]))]
		data = [ SCONFIG[1][d][1] for d in range(n_depth)]
		positions = [ i+s+1 for i in base_pos ]
		col = SCONFIG[0]
		bp = plt.boxplot(data, positions=positions, sym='', widths=BOX_WIDTH)
		#for i in bp['boxes']:
		#	i.set_facecolor(col)
		plt.setp(bp['boxes'], color = col)
		plt.setp(bp['whiskers'], color = col)
		plt.setp(bp['caps'], color = col)
		plt.setp(bp['medians'], color = "black")
	plt.xticks([ i + (float(n_strains+1)/2.0) for i in base_pos ], xlabels)
	plt.xlim(0,(n_strains+1)*n_depth)
	#fig.set_size_inches(14,8)
	plt.ylabel("Time to execute DeletionID\n(in seconds)")
	plt.xlabel("Dataset size\n(# of paired-end reads sampled)")
	plt.title(args.title)
	#plt.show()
	out_png_fn = args.title.replace(" ","_")+".png"
	plt.savefig(out_png_fn)
	print(out_png_fn)
