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

roman2arabic = {"chrI":"chr1","chrII":"chr2","chrIII":"chr3","chrIV":"chr4","chrV":"chr5",
			"chrVI":"chr6","chrVII":"chr7","chrVIII":"chr8","chrIX":"chr9","chrX":"chr10",
			"chrXI":"chr11","chrXII":"chr12","chrXIII":"chr13","chrXIV":"chr14","chrXV":"chr15",
			"chrXVI":"chr16",}

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-t','--title', metavar='figure_title', dest='title', required=True, help='')
	parser.add_argument('-c','--config', metavar='config_fn', dest='config_fn', required=True, help='the config file for a grid-organized subplot')
	parser.add_argument('-g','--features-gff', metavar='features_gff', dest='features_gff', required=True, help='the featuer GFF file from SGD to get the gene coordinates')
	parser.add_argument('-d','--header', dest='header', default=False, required=False, help='skip first line as column header')

	args = parser.parse_args()
	return(args)

def get_bedgraph_info(bedgraph_file, locus_coord, flanking):
	flanking_range = (roman2arabic[locus_coord[0]], locus_coord[1]-flanking, locus_coord[2]+flanking)
	x_vector = list(range(flanking_range[1],flanking_range[2]))
	y_vector = [0] * len(x_vector)
	reader = open(bedgraph_file,'r')
	for line in reader:
		tokens = line.strip().split('\t')
		# Skip if chromosome doesn't match
		if(tokens[0]!=flanking_range[0]):
			continue
			# Skip if interval before interval of interest
		elif(int(tokens[1])<flanking_range[1] and int(tokens[2])<flanking_range[1]):
			continue
			# Skip if interval after interval of interest
		elif(int(tokens[1])>flanking_range[2] and int(tokens[2])>flanking_range[2]):
			continue
		value = int(tokens[3])
		for local_x in range(int(tokens[1]),int(tokens[2])):
			if(local_x in x_vector):
				y_vector[x_vector.index(local_x)] = value
	reader.close()
	plot_info = {"X":x_vector,"Y":y_vector,"chrom":locus_coord[0]}
	return(plot_info)

def parse_configs(configs_fn):
	subplot_configs = {}
	reader = open(configs_fn,'r')
	for line in reader:
		tokens = line.strip().split('\t')
		if(tokens[0]=="FLANKING"):
			subplot_configs.update({tokens[0]:int(tokens[1])})
			continue
		elif(tokens[0]=="S_MAX_Y"):
			subplot_configs.update({tokens[0]:[ int(i) for i in tokens[1:]]})
			continue
		subplot_configs.update({tokens[0]:tokens[1:]})
	reader.close()
	# Count subplot dimensions
	subplot_configs.update({"N_SAMPLES":len(subplot_configs["SAMPLES"])})
	subplot_configs.update({"N_LOCI":len(subplot_configs["LOCI"])})
	# Validate configs
	for key in ["S_LABEL","S_COLOR","S_MAX_Y"]:
		if(len(subplot_configs[key])!=subplot_configs["N_SAMPLES"]):
			sys.stderr.write("Mismatch in number of samples with %i field. Exiting...\n" % (key))
			quit()
	return(subplot_configs)

def parse_gff(gff_fn, loci_list):
	locus2coord = {}
	reader = open(gff_fn,'r')
	for line in reader:
		if(line.find("#")==0):
			continue
		if(line.find(">")==0):
			break
		tokens = line.strip().split('\t')
		gene_name = ""
		for feature in tokens[8].split(';'):
			if(feature.find("gene=")!=0):
				continue
			gene_name = feature.split('=')[1]
			break
		if(gene_name in loci_list):
			locus2coord.update({gene_name:(tokens[0],int(tokens[3])-1,int(tokens[4]))})
	reader.close()
	return(locus2coord)

if __name__ == "__main__":
	'''Plot scatter'''
	args = getParams()

	CONFIGS = parse_configs(args.config_fn)
	LOCUS2COORD = parse_gff(args.features_gff, CONFIGS["LOCI"])

	fig, asx = plt.subplots(CONFIGS["N_SAMPLES"],CONFIGS["N_LOCI"])
	fig.suptitle(args.title)
	plt.tight_layout()
	for s in range(CONFIGS["N_SAMPLES"]):
		bedgraph_fn = "results/BedGraphs/%s.raw.bedgraph" % CONFIGS["SAMPLES"][s]
		for l in range(CONFIGS["N_LOCI"]):
			locus = CONFIGS["LOCI"][l]
			sys.stderr.write("Processing sample %s by locus %s...\n" % (bedgraph_fn, locus))
			data = get_bedgraph_info(bedgraph_fn, LOCUS2COORD[locus], CONFIGS["FLANKING"])
			# Plot data
			asx[s,l].fill_between(data["X"], data["Y"], color=CONFIGS["S_COLOR"][s])
			asx[s,l].set_ylim(bottom=0,top=CONFIGS["S_MAX_Y"][s])
			asx[s,l].label_outer()
			x0 = data["X"][0]
			xend = data["X"][-1]
			xstart = data["X"][0] + CONFIGS["FLANKING"]
			xstop = data["X"][-1] - CONFIGS["FLANKING"]
			asx[s,l].set_xlim([x0,xend])
			plt.sca(asx[s,l])
			plt.xticks([x0,xstart,xstop,xend],["-200","start","stop","+200"])
		# Label Samples
		asx[s,0].set_ylabel(CONFIGS["S_LABEL"][s])

	for l in range(CONFIGS["N_LOCI"]):
		coord = LOCUS2COORD[CONFIGS["LOCI"][l]]
		asx[CONFIGS["N_SAMPLES"]-1,l].set_xlabel("%s:%i-%i" % (coord[0],coord[1],coord[2]))
	fig.set_size_inches(14,8)
	#plt.show()
	out_pic_fn = args.title.replace(" ","_")+".svg"
	plt.savefig(out_pic_fn)
	print(out_pic_fn)
