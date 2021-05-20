from os import listdir
from os.path import isfile, join
import sys
import re
import argparse

# Python 3.6+
# relies on dict insertion order
roman2arabic = {"chrI":"chr1","chrII":"chr2","chrIII":"chr3","chrIV":"chr4","chrV":"chr5",
			"chrVI":"chr6","chrVII":"chr7","chrVIII":"chr8","chrIX":"chr9","chrX":"chr10",
			"chrXI":"chr11","chrXII":"chr12","chrXIII":"chr13","chrXIV":"chr14","chrXV":"chr15",
			"chrXVI":"chr16",}



def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Use pileup information to get a heatmap of each sample\'s coverage at the expected KO site.')

	parser.add_argument('-s','--summary_fn', metavar='summary_fn', dest='summary_fn', required=True, help='the output_summary file failed entries')
	parser.add_argument('-i','--input-dir', metavar='input_dir', dest='input_dir', required=True, help='the directory where all the BedGraph pileup files are saved')
	parser.add_argument('-g','--features-gff', metavar='features_gff', dest='features_gff', required=True, help='the featuer GFF file from SGD to get the gene coordinates')
	parser.add_argument('-w','--window', metavar='win_size', dest='window', default=6000, type=int, help='the window size to center around feature range')

	args = parser.parse_args()
	return(args)

#chr2   334391  336812  REB1    .   -
def parse_bedfile(bed_fn):
	orf2coord = {}
	reader = open(bed_fn,'r')
	for line in reader:
		tokens = line.strip().split('\t')
		orf2coord.update({tokens[3]:(tokens[0],int(tokens[1]),int(tokens[2]))})
	reader.close()
	return(orf2coord)

def expand_coord(bed_coord, window):
	midpoint = bed_coord[1] + (bed_coord[2] - bed_coord[1])//2
	flank = window//2
	return(bed_coord[0],midpoint-flank,midpoint+flank)

# chr1	0	1	14
# chr1	1	2	20
# chr1	2	3	25
# chr1	3	5	27
def parse_bedgraph(bg_fn, bed_coord, window=2000):
	window_coord = expand_coord(bed_coord,window)
	pileup = ["NaN"] * (window_coord[2]-window_coord[1])
	reader = open(bg_fn, 'r')
	for line in reader.readlines():
		if(line.find('#')==0):
			continue
		tokens = line.strip().split('\t')
		if(tokens[0]!=window_coord[0]):
			continue
			# Skip if interval before interval of interest
		elif(int(tokens[1])<window_coord[1] and int(tokens[2])<window_coord[1]):
			continue
			# Skip if interval after interval of interest
		elif(int(tokens[1])>window_coord[2] and int(tokens[2])>window_coord[2]):
			continue
		value = float(tokens[3])
		for local_x in range(int(tokens[1]),int(tokens[2])):
			if(local_x >= window_coord[1] and local_x<window_coord[2]):
				pileup[local_x-window_coord[1]] = value
	reader.close()
	return(pileup)


def parse_gff(gff_fn):
	locus2coord = {}
	reader = open(gff_fn,'r')
	for line in reader:
		if(line.find("#")==0):
			continue
		if(line.find(">")==0):
			break
		tokens = line.strip().split('\t')
		if(tokens[2] in ["mRNA","CDS"]):
			continue
		gene_name = ""
		for feature in tokens[8].split(';'):
			if(feature.find("ID=")!=0):
				continue
			gene_name = feature.split('=')[1]
			break
		locus2coord.update({gene_name:(roman2arabic.get(tokens[0],"chrZ"),int(tokens[3])-1,int(tokens[4]))})
	reader.close()
	return(locus2coord)


#STATUS	feature_type	NOTES	KO_SCORE	SYS	STD	TableS1_Deletion	TableS1_replicate_id	ERS_accession	n_hits	hit_list	hit_scores	LEU2_SCORE	URA3_SCORE	experiment_accession	run_accession	submission_accession	nominal_length	read_count	base_count	first_public	nominal_sdev
#PASS	ORF-Uncharacterized		-6.109952403138639	YAL064C-A	TDA8	Del1_TDA8	SD0863b	ERS838232	3	LEU2|URA3|TDA8	ND|ND|-6.109952403138639	ND	ND	ERX1406336	ERR1334744	ERA587837	484	8807338	1329908038	2016-03-22	81
#PASS	ORF-Uncharacterized		-5.807910468072448	YAL064C-A	TDA8	Del1_TDA8	SD0863b2	ERS838233	3	LEU2|URA3|TDA8	ND|ND|-5.807910468072448	ND	ND	ERX1406337	ERR1334745	ERA587837	484	8996386	1358454286	2016-03-22	81
#FAIL	ORF-Verified		-	YBL091C-A	SCS22	Del2_SCS22	SD0864b	ERS838234	2	LEU2|URA3	ND|ND	ND	ND	ERX1406338	ERR1334746	ERA587837	484	8710346	1315262246	2016-03-22	81
#FAIL	ORF-Verified		-	YBL091C-A	SCS22	Del2_SCS22	SD0864b2	ERS838235	2	LEU2|URA3	ND|ND	ND	ND	ERX1406339	ERR1334747	ERA587837	484	8579514	1295506614	2016-03-22	81
if __name__ == "__main__":
	'''Collect metadata and DeletionID results to get detection stats on the YKOC data'''

	hardcode_name_remap = {
		"YCR061W":"TVS1",
		"YCR100C":"EMA35",
		"YFR045W":"MRX20",
		"YIR035C":"NRE1",
		"YNR062C":"PUL3",
		"YNR063W":"PUL4",
		"YER156C":"MYG1",
		"YMR087W":"PDL32",
		"YLR050C":"EMA19",
		"YMR279C":"ATR2",
		"YMR102C":"LAF1",
		"YMR111C":"EUC1",
		"YMR130W":"DPI35",
		"YJR039W":"MLO127",
		"YJR061W":"MNN14",
		"YGR053C":"MCO32",
		"YKR023W":"RQT4",
		"PET10":"PLN1"
	}

	# Get params
	args = getParams()
	WINDOW = args.window
	orf2bed = parse_gff(args.features_gff)
	ers2pileup = {}

	# Parse metadata
	reader = open(args.summary_fn, 'r')#, encoding='utf-8')
	for mline in reader:
		mtokens = mline.strip().split('\t')
		ERS = mtokens[9]
		#ORF = hardcode_name_remap.get(mtokens[5],mtokens[5])
		SYS = mtokens[5]
		STD = mtokens[6]
		COORD = orf2bed.get(SYS,("chrZ",0,1))

		sys.stderr.write("Parsing for sample %s\n" % ERS )

		# Build BedGraph filename
		# Pileup BedGraph in CDT interval, then normalize
		bedgraph_fn = join(args.input_dir,"%s.bedgraph" % ERS)
		if(isfile(bedgraph_fn)):
			pileup = parse_bedgraph(bedgraph_fn, COORD, WINDOW)
		else:
			sys.stderr.write("%s BedGraph pileup does not exist\n" % bedgraph_fn)
			continue

		# Build CDT filename
		cdt_row_fn = join(args.input_dir,"%s_%s_%ibp.cdt" % (ERS,STD,WINDOW))
		writer = open(cdt_row_fn,'w')
		# Write CDT header
		writer.write("\t".join([ "YORF", "NAME", "LENGTH" ]) + "\t" + \
					"\t".join([ str(i) for i in range(WINDOW)]) + "\n")
		# Write CDT row vals
		writer.write("\t".join([ ERS, STD, str(COORD[2]-COORD[1]) ]) + "\t" + \
					"\t".join([ str(i) for i in pileup]) + "\n")
		writer.close()

		# (ERS,orf_len):[pileup CDT row]
		# cdt_line = "%s\t%s\t%s" % (ERS,ORF,"\t".join([str(i) for i in pileup]))
		# ers2pileup.update({(ERS,COORD[2]-COORD[1]):cdt_line})

		sys.stderr.write("%s pileup complete\n" % ERS)

	reader.close()

	# # Write CDT header
	# sys.stdout.write("\t".join([ "YORF", "NAME"]) + "\t" + \
	# 					"\t".join([ str(i) for i in range(WINDOW)]) + "\n")
	# # Write Output by gene length
	# COORD_KEYS = sorted(ers2pileup, key=lambda x:x[1])
	# sys.stdout.write("\n".join([ers2pileup[KEY] for KEY in COORD_KEYS]))
