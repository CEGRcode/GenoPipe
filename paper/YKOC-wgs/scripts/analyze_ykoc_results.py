from os import listdir
from os.path import isfile, join
import sys
import re
import argparse

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Parse several metadata files and DeletionID output to check detection rates of the DeletionID module on the YKOC data from Puddu et al.')
	parser.add_argument('-m','--metadata', metavar='metadata_fn', dest='metadata_fn', required=True, help='the metadata file downloaded of Table S1 from the paper that maps deletion background to the ERS accession')
	parser.add_argument('-i','--input-dir', metavar='input_dir', dest='input_dir', required=True, help='the directory where all the DeletionID output files were saved to (*_deletion.tab)')
	
	parser.add_argument('-e','--ebi-metadata', metavar='ebi_metadata_fn', dest='ebi_metadata_fn', default="", required=False, help='the file with EBI metadata includin read counts')
	parser.add_argument('-t','--verifyKO-info', metavar='table_s6_file', dest='tables6_fn', required=True, help='the file with information on verified knockouts from Puddu et al')
	
	parser.add_argument('-n','--name-map', metavar='names_fn', dest='names_fn', default="", required=False, help='the GFF file with mappings of SGD standard and systematic names to aliases')
	parser.add_argument('-d','--deleted-map', metavar='deleted_features_file', dest='deleted_fn', required=True, help='the file with information on deleted and merged SGD annotations')

	args = parser.parse_args()
	return(args)

#YAR043C ORF|Deleted 1   193248  192928  C   S000120182              Deleted ORF; does not encode a protein; included in the original annotation of Chromosome I but later withdrawn Deleted ORF, does not encode a protein; this putative ORF was included in the original annotation of Chromosome I but was later withdrawn   1999-07-17
#YAR044W ORF|Merged  1   193600  196179  W   S000000082  L000001316  YAR042W S000000081  Merged open reading frame; does not encode a discrete protein; YAR044W was originally annotated as an independent ORF, but as a result of a sequence change, it was merged with an adjacent ORF into a single reading frame, designated YAR042W     2003-09-27
#YAR052C ORF|Deleted 1   207502  207122  C   S000120184              Deleted ORF; does not encode a protein; included in the original annotation of Chromosome I but later withdrawn Deleted ORF, does not encode a protein; this putative ORF was included in the original annotation of Chromosome I but was later withdrawn   1999-07-17
#YAR062W pseudogene|Merged   1   218549  219145  W   S000000088      YAR061W S000000087  Pseudogenic fragment with similarity to flocculins; YAR062W has been merged into YAR061W; this region has a paralog, YHR213W, that arose from a segmental duplication   As part of SGD's genome annotation revision R64.2, the two pseudogenes YAR061W + YAR062W have been combined into a single pseudogene, keeping the name YAR061W. This region is homologous to parts of FLO1/YAR050W and FLO9/YAL063C, but contains stop codons at several positions. 2014-11-18
def get_deleted_map(deleted_fn):
	'''Collect names of deleted or merged genes and map to their deleted or merged status'''
	deletion_map = {"YAR037W":("Undocumented",""),"YAR040C":("Undocumented","")}
	reader = open(deleted_fn,'r')
	for line in reader:
		tokens = line.strip().split('\t')
		key = tokens[0]
		if(tokens[1].split('|')[0] not in ["pseudogene","ORF"]):
			continue
		#print("Skipping type=%s" % tokens[1].split('|')[0])
		status = tokens[1].split('|')[1]
		alias = tokens[8]
		#print("status: %s, alias: %s" % (status,alias))
		deletion_map[key] = (status,alias)
	reader.close()
	return(deletion_map)

#chrI    SGD gene    335 649 .   +   .   ID=YAL069W;Name=YAL069W;Ontology_term=GO:0003674,GO:0005575,GO:0008150,SO:0000704;Note=Dubious%20open%20reading%20frame%3B%20unlikely%20to%20encode%20a%20functional%20protein%2C%20based%20on%20available%20experimental%20and%20comparative%20sequence%20data;display=Dubious%20open%20reading%20frame;dbxref=SGD:S000002143;orf_classification=Dubious;curie=SGD:S000002143
#chrI    SGD gene    538 792 .   +   .   ID=YAL068W-A;Name=YAL068W-A;Ontology_term=GO:0003674,GO:0005575,GO:0008150,SO:0000704;Note=Dubious%20open%20reading%20frame%3B%20unlikely%20to%20encode%20a%20functional%20protein%2C%20based%20on%20available%20experimental%20and%20comparative%20sequence%20data%3B%20identified%20by%20gene-trapping%2C%20microarray-based%20expression%20analysis%2C%20and%20genome-wide%20homology%20searching;display=Dubious%20open%20reading%20frame;dbxref=SGD:S000028594;orf_classification=Dubious;curie=SGD:S000028594
#chrI    SGD gene    1807    2169    .   -   .   ID=YAL068C;Name=YAL068C;gene=PAU8;Alias=PAU8,seripauperin%20PAU8;Ontology_term=GO:0003674,GO:0005575,GO:0030437,GO:0045944,SO:0000704;Note=Protein%20of%20unknown%20function%3B%20member%20of%20the%20seripauperin%20multigene%20family%20encoded%20mainly%20in%20subtelomeric%20regions;display=Protein%20of%20unknown%20function;dbxref=SGD:S000002142;orf_classification=Verified;curie=SGD:S000002142
#chrI    SGD gene    2480    2707    .   +   .   ID=YAL067W-A;Name=YAL067W-A;Ontology_term=GO:0003674,GO:0005575,GO:0008150,SO:0000704;Note=Putative%20protein%20of%20unknown%20function%3B%20identified%20by%20gene-trapping%2C%20microarray-based%20expression%20analysis%2C%20and%20genome-wide%20homology%20searching;display=Putative%20protein%20of%20unknown%20function;dbxref=SGD:S000028593;orf_classification=Uncharacterized;curie=SGD:S000028593
#chrI    SGD gene    7235    9016    .   -   .   ID=YAL067C;Name=YAL067C;gene=SEO1;Alias=SEO1,putative%20permease%20SEO1;Ontology_term=GO:0015124,GO:0016020,GO:0055085,GO:0055085,GO:0071944,SO:0000704;Note=Putative%20permease%3B%20member%20of%20the%20allantoate%20transporter%20subfamily%20of%20the%20major%20facilitator%20superfamily%3B%20mutation%20confers%20resistance%20to%20ethionine%20sulfoxide;display=Putative%20permease;dbxref=SGD:S000000062;orf_classification=Verified;curie=SGD:S000000062
#chrI    SGD gene    10091   10399   .   +   .   ID=YAL066W;Name=YAL066W;Ontology_term=GO:0003674,GO:0005575,GO:0008150,SO:0000704;Note=Dubious%20open%20reading%20frame%3B%20unlikely%20to%20encode%20a%20functional%20protein%2C%20based%20on%20available%20experimental%20and%20comparative%20sequence%20data;display=Dubious%20open%20reading%20frame;dbxref=SGD:S000000061;orf_classification=Dubious;curie=SGD:S000000061
#chrI    SGD gene    11565   11951   .   -   .   ID=YAL065C;Name=YAL065C;Ontology_term=GO:0003674,GO:0005575,GO:0008150,SO:0000704;Note=Putative%20protein%20of%20unknown%20function%3B%20shows%20sequence%20similarity%20to%20FLO1%20and%20other%20flocculins;display=Putative%20protein%20of%20unknown%20function;dbxref=SGD:S000001817;orf_classification=Uncharacterized;curie=SGD:S000001817
#chrI    SGD gene    12046   12426   .   +   .   ID=YAL064W-B;Name=YAL064W-B;Ontology_term=GO:0003674,GO:0005783,GO:0008150,SO:0000704;Note=Fungal-specific%20protein%20of%20unknown%20function;display=Fungal-specific%20protein%20of%20unknown%20function;dbxref=SGD:S000002141;orf_classification=Uncharacterized;curie=SGD:S000002141
#chrI    SGD gene    13363   13743   .   -   .   ID=YAL064C-A;Name=YAL064C-A;gene=TDA8;Alias=TDA8,YAL065C-A;Ontology_term=GO:0003674,GO:0005575,GO:0008150,SO:0000704;Note=Putative%20protein%20of%20unknown%20function%3B%20null%20mutant%20is%20sensitive%20to%20expression%20of%20the%20top1-T722A%20allele%3B%20not%20an%20essential%20gene;display=Putative%20protein%20of%20unknown%20function;dbxref=SGD:S000002140;orf_classification=Uncharacterized;curie=SGD:S000002140
def get_orf_info(orf_info_fn):
	'''Create a mapping of alias names to systematic gene names'''
	alias2sys = {"DUR1,2":"YBR208C","ADE5,7":"YGL234W","ARG5,6":"YER069W","MFALPHA1":"YPL187W","MFALPHA2":"YGL089C","HUF1":"YOR300W"}
	#"ENV10" }
	sys2info = {}
	reader = open(orf_info_fn,'r')
	for line in reader:
		if(line.find("#")==0):
			continue
		elif(line.find(">")==0):
			break
		tokens = line.strip().split('\t')
		feature_type = tokens[2]
		#print(feature_type)
		if(feature_type not in ["gene","ORF","pseudogene","blocked_reading_frame","transposable_element_gene","mating_type_region"]):
			continue
		sys = "SYSTEMATIC"
		std = "STANDARD"
		aliases = []
		for attribute in tokens[8].split(';'):
			if(attribute.find("ID=")==0):
				sys = attribute.split("=")[1]
			elif(attribute.find("gene=")==0):
				std = attribute.split("=")[1]
			elif(attribute.find("Alias=")==0):
				aliases = attribute.split("=")[1].split(",")
			elif(attribute.find("orf_classification=")==0):
				feature_type += "-%s" % (attribute.split("=")[1])
		if(std=="STANDARD"):
			std = sys
		alias2sys.update({sys:sys,std:sys})
		for a in aliases:
			if(re.fullmatch("[A-Z]{3}[0-9]+",a) or re.fullmatch("Y[A-P][RL][0-9]{3}[WC](-[ABCD])+",a)):
				alias2sys.setdefault(a, sys)
		sys2info.update({sys:(std,feature_type)})
	reader.close()
	return(alias2sys,sys2info)


#LEU2    No Data Detected
#URA3    No Data Detected
#CHO2    -3.864576691071892
def get_id_hits_and_scores(id_file):
	'''Parse id file for target info & epitope read count'''
	id_genes = {}
	id_reader = open(id_file,'r')
	for line in id_reader:
		tokens = line.strip().split('\t')
		if(tokens[1]=="No Data Detected"):
			tokens[1] = "ND"
		id_genes.update({tokens[0]:tokens[1]})
	id_reader.close()
	return(id_genes)

#sample_accession	secondary_sample_accession	experiment_accession	run_accession	submission_accession	nominal_length	read_count	base_count	first_public	fastq_bytes	fastq_md5	fastq_ftp	fastq_aspera	submitted_bytes	submitted_md5	submitted_ftp	submitted_aspera	nominal_sdev
#SAMEA3531083	ERS838232	ERX1406336	ERR1334744	ERA587837	484	8807338	1329908038	2016-03-22	203589099;292667412	88995985b66397b5fb0743ccea16b167;9ef867bd3a25b7343454f18faece0179	ftp.sra.ebi.ac.uk/vol1/fastq/ERR133/004/ERR1334744/ERR1334744_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR133/004/ERR1334744/ERR1334744_2.fastq.gz	fasp.sra.ebi.ac.uk:/vol1/fastq/ERR133/004/ERR1334744/ERR1334744_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/ERR133/004/ERR1334744/ERR1334744_2.fastq.gz	333265131	8e1a1414960c51d41c298635c773a4dbftp.sra.ebi.ac.uk/vol1/run/ERR133/ERR1334744/18596_1#1.cram	fasp.sra.ebi.ac.uk:/vol1/run/ERR133/ERR1334744/18596_1#1.cram	81
#SAMEA3531084	ERS838233	ERX1406337	ERR1334745	ERA587837	484	8996386	1358454286	2016-03-22	207775722;293057592	85f615bc0598404b0da389bcd26af18a;d4fcb005f1e906221bde0557cb87db6c	ftp.sra.ebi.ac.uk/vol1/fastq/ERR133/005/ERR1334745/ERR1334745_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR133/005/ERR1334745/ERR1334745_2.fastq.gz	fasp.sra.ebi.ac.uk:/vol1/fastq/ERR133/005/ERR1334745/ERR1334745_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/ERR133/005/ERR1334745/ERR1334745_2.fastq.gz	336752949	7965c11ea7e13605bca5257b14bd35c7ftp.sra.ebi.ac.uk/vol1/run/ERR133/ERR1334745/18596_1#2.cram	fasp.sra.ebi.ac.uk:/vol1/run/ERR133/ERR1334745/18596_1#2.cram	81
def get_ebi_metadata(ebi_metadata_fn):
	'''Parse EBI metadata to report this info keyed on the ERS code'''
	ebi_info = {}
	reader = open(ebi_metadata_fn)
	for line in reader:
		tokens = line.strip().split('\t')
		subset_tokens = tokens[2:9]
		subset_tokens.append(tokens[-1])
		ebi_info.update({tokens[1]:'\t'.join(subset_tokens)})
	reader.close()
	return(ebi_info)

###DNTAG_freqs	Percentage of times DNTAG_seqs have been detected between D1 and D2 sequences (comma separated list)							
###DNTAG_notes	Notes (NA indicates detected DNTAG matches expected DNTAG)							
##DelName	SDname	Confirmed_deletion	UPTAG_seqs	UPTAG_freqs	UPTAG_notes	DNTAG_seqs	DNTAG_freqs	DNTAG_notes
#Del5_UBX7	SD0867b	YBR273C	ATTTGTGTCCCTGCGCCTG	92.59	Tag is different from reference expected tag	GGCGATATCCAGCCTAG	97.44	Tag is different from reference expected tag
#Del5_UBX7	SD0867b2	YBR273C	ATTTGTGTCCCTGCGCCTG	97.14	Tag is different from reference expected tag	GGCGATATCCAGCCTAG	97.67	Tag is different from reference expected tag
#Del182_GPT2	SD1044b	YKR067W	ATACGCTGTCGGAATACCG	93.75	Tag is different from reference expected tag	ATTCTTGCAGACCAGAGCG	92	Tag is different from reference expected tag
#Del182_GPT2	SD1044b2	YKR067W	ATACGCTGTCGGAATACCG	100	Tag is different from reference expected tag	ATTCTTGCAGACCAGAGCG	93.94	Tag is different from reference expected tag
#Del430_MET6	SD1292b	YER091C	AGTGCGGCCATCGTCATTAT	83.33	Tag is different from reference expected tag	AGCTCATCCAACGCTTATA	96.97	Tag is different from reference expected tag
#Del430_MET6	SD1292b2	YER091C	AGTGCGGCCATCGTCATTAT	93.75	Tag is different from reference expected tag	AGCTCATCCAACGCTTATA	100	Tag is different from reference expected tag
def get_verifyKO_info(puddu_table_s6_fn):
	'''Parse Puddu et al 2019 Table S6 for verification info on each knockout keyed on SD id (id system used by the paper)'''
	ts6_info = {}
	reader = open(puddu_table_s6_fn)
	for line in reader:
		tokens = line.strip().split('\t')
		if(line.find("##")==0):
			continue
		elif(line.find("#")==0):
			ts6_info.update({tokens[1]:"\t".join(tokens[2:])})
			continue
		ts6_info.update({tokens[1]:"\t".join(tokens[2:])})
	reader.close()
	return(ts6_info)

###Supplementary Table 1. Accession numbers of YKOC colonies sequenced.
###LEGEND
###Deletion_Strain       Unique ID identifying the KO strain from which colonies were derived
###SD_code       Unique ID identifying a single colony
###ERS accession number  Accession code to download data from the European Nucleotide Archive
##Deletion_Strain        SD_code ERS accession number
#Del1_TDA8       SD0863b ERS838232
#Del1_TDA8       SD0863b2        ERS838233
#Del2_SCS22      SD0864b ERS838234
#Del2_SCS22      SD0864b2        ERS838235
if __name__ == "__main__":
	'''Collect metadata and DeletionID results to get detection stats on the YKOC data by iterating through Puddu et al Table S1'''
	# Get params
	args = getParams()
	# Get EBI info map
	ers2ebi_info = get_ebi_metadata(args.ebi_metadata_fn)
	# Get Table S6 info
	sd2verifyKO_info = get_verifyKO_info(args.tables6_fn)
	# Initialize deletion map
	orf2deleted = get_deleted_map(args.deleted_fn)
	# Initialize alias map
	name2sys,sys2info = get_orf_info(args.names_fn)
	# Initialize tracking variables
	found = 0
	total = 0
	# Write header
	# STATUS, NOTES, KO_SCORE, SYS, STD,
	sys.stdout.write("\t".join([ "STATUS", "feature_type", "TableS6_verified", "NOTES", "KO_SCORE", "SYS", "STD", "TableS1_Deletion", "SD_code", "ERS_accession", \
			"n_hits", "hit_list", "hit_scores", \
			"LEU2_SCORE","URA3_SCORE", \
			ers2ebi_info["secondary_sample_accession"], sd2verifyKO_info["SDname"]]) + "\n")
	
	# Parse metadata
	reader = open(args.metadata_fn, 'r', encoding='utf-8')
	for mline in reader:
		if(mline.find('#')==0):
			continue
		mtokens = mline.strip().split('\t')
		
		# Status ultimately set to "OMIT" "PASS" and "FAIL"
		status = "UNKNOWN"
		notes = []
		ko_score = "-"
		id_genes = {}
		
		# Pull relevant info from metadata tokens
		sd_name = mtokens[1]
		ers = mtokens[2]
		ebi_info = ers2ebi_info.get(ers,'\t'.join(["NoInfo"]*8))
		ts6_info = sd2verifyKO_info.get(sd_name,"\t".join([""]*7))
		ts6_verified = "NotVerified"
		if(sd_name in sd2verifyKO_info.keys()):
			ts6_verified = "Verified"
		read_count = -1
		if(ebi_info.find("NoInfo")!=0):
			read_count = int(ebi_info.split('\t')[4])
		ko_gene = mtokens[0].split('_')[1]
		SYS = name2sys.get(ko_gene, "DNE")
		STD = sys2info.get(SYS, ("-","TypeNotAvailable"))[0]
		TYPE = sys2info.get(SYS, ("-","TypeNotAvailable"))[1]
		
		# Parse id file (if exists) for target info & epitope read count
		id_file = join(args.input_dir,"%s_deletion.tab" % ers)
		if(isfile(id_file)):
			id_gene2score = get_id_hits_and_scores(id_file)
		else:
			notes.append("ID-file-doees-not-exist")
		
		if(read_count < 500000):
			notes.append("less-than-500K")
		elif(read_count < 1000000):
			notes.append("less-than-1M")
		elif(read_count < 2000000):
			notes.append("less-than-1M")
		
		# Check if Deleted/Merged
		if(ko_gene in orf2deleted.keys()):
			TYPE = orf2deleted[ko_gene][0]
			status = "OMIT"
			notes.append('-'.join(orf2deleted[ko_gene]))
		# Check if feature type excluded from DeletionID db
		elif(TYPE == "pseudogene" or TYPE == "transposable_element_gene"):
			status = "OMIT"
		# Check if Wildtype control
		elif(ko_gene in ["WT-1","WT-2","WT-3","WT-4"]):
			status = "CONTROL"
			TYPE = "-"
			notes.append("WT-control")
		# Check if ko in list of hits
		elif( SYS in [name2sys.get(k,"-") for k in id_gene2score.keys()]):
			status = "PASS"
			found += 1
			ko_score = id_gene2score.get(STD,"ScoreMapFail")
		else:
			status = "FAIL"
		total += 1
		# Output:
		sys.stdout.write( "\t".join([ status, TYPE, ts6_verified, "|".join(notes), ko_score, SYS, STD, "\t".join(mtokens), \
				str(len(id_gene2score.keys())), "|".join(id_gene2score.keys()), "|".join(id_gene2score.values()), \
				id_gene2score.get("LEU2","-"), id_gene2score.get("URA3","-"), \
				ebi_info, ts6_info ]) + "\n")
	reader.close()
	# Output Summary
	sys.stderr.write('Output Summary:\n')
	sys.stderr.write('%i knock-outs found of %i samples\n' % (found, total))
