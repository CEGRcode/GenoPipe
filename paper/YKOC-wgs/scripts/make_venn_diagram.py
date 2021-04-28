import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Generates Figure3C, a proportional Venn diagram of DeletionID\'s successful KO identification vs the Puddu et al method')

	parser.add_argument('-t','--title', metavar='figure_title', dest='title', required=True, help='')
	parser.add_argument('-s','--summary_fn', metavar='summary_fn', dest='summary_fn', required=True, help='the ykoc_output_summary.txt file')

	args = parser.parse_args()
	return(args)

#STATUS	feature_type	TableS6_verified	NOTES	KO_SCORE	SYS	STD	TableS1_Deletion	SD_code	ERS_accession	n_hits	hit_list	hit_scores	LEU2_SCORE	URA3_SCORE	experiment_accession	run_accession	submission_accession	nominal_length	read_count	base_count	first_public	nominal_sdev	Confirmed_deletion	UPTAG_seqs	UPTAG_freqs	UPTAG_notes	DNTAG_seqs	DNTAG_freqs	DNTAG_notes
#PASS	gene-Uncharacterized	Verified		-5.679247695	YAL064C-A	TDA8	Del1_TDA8	SD0863b	ERS838232	3	LEU2|URA3|TDA8	ND|ND|-5.6792476952599165	ND	ND	ERX1406336	ERR1334744	ERA587837	484	8807338	1329908038	3/22/16	81	YAL064C-A	CCTGATGAGTAATGCGATAG	100	NA	ACACGATGATAAGCCAGCGC	89.47	NA
#PASS	gene-Uncharacterized	Verified		-5.681115917	YAL064C-A	TDA8	Del1_TDA8	SD0863b2	ERS838233	3	LEU2|URA3|TDA8	ND|ND|-5.681115917493542	ND	ND	ERX1406337	ERR1334745	ERA587837	484	8996386	1358454286	3/22/16	81	YAL064C-A	CCTGATGAGTAATGCGATAG	92.5	NA	ACACGATGATAAGCCAGCGC	100	NA
#FAIL	gene-Verified	Verified		-	YBL091C-A	SCS22	Del2_SCS22	SD0864b	ERS838234	2	LEU2|URA3	ND|ND	ND	ND	ERX1406338	ERR1334746	ERA587837	484	8710346	1315262246	3/22/16	81	YBL091C-A	TCATAGTACGCATGTCGAGG	94.59	NA	ATTGTAACGACGCGCTTGC	100	Tag is different from reference expected tag
#FAIL	gene-Verified	Verified		-	YBL091C-A	SCS22	Del2_SCS22	SD0864b2	ERS838235	2	LEU2|URA3	ND|ND	ND	ND	ERX1406339	ERR1334747	ERA587837	484	8579514	1295506614	3/22/16	81	YBL091C-A	TCATAGTACGCATGTCGAGG	97.44	NA	ATTGTAACGACGCGCTTGC	90	Tag is different from reference expected tag
if __name__ == "__main__":
	'''Parse metadata to get a confirmed KO venn diagram'''
	# Get params
	args = getParams()

	# Parse metadata
	p_group = []
	d_group = []
	reader = open(args.summary_fn, 'r')#, encoding='utf-8')
	for mline in reader:
		mtokens = mline.strip().split('\t')
		DELID_KO = mtokens[0]
		PUDDU_KO = mtokens[2]
		ERS = mtokens[9]
		# Skip header, controls, and omitted KOs
		if(DELID_KO in ["OMIT","CONTROL","STATUS"]):
			continue
		# Tally Puddu et al verified
		if(PUDDU_KO == "Verified"):
			p_group.append(ERS)
		# Tally DeletionID verified
		if(DELID_KO=="PASS"):
			d_group.append(ERS)
	reader.close()

	venn2([set(p_group),set(d_group)], set_labels = ("KO confirmed by\nPuddu et al 2019","KO confirmed by\nDeletionID"), \
			set_colors = ("red","blue"))
	plt.title(args.title)
	#plt.show()
	out_pic_fn = args.title.replace(" ","_")+".svg"
	plt.savefig(out_pic_fn)
	print(out_pic_fn)
