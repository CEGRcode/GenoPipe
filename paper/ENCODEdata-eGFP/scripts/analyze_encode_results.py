from os import listdir
from os.path import isfile, join
import sys
import argparse

# Python 3 needed for encoding feature for UTF-8
# (ENCODE uses some capital delta chars in summary descriptions of GeneticModifications)

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Parse metadata file and GenoPipe output to check detection rates of the GenoPipe tool.')
	parser.add_argument('-m','--metadata', metavar='metadata_fn', required=True, help='the metadata file downloaded with ENCODE dataset that includes info like PE/SE, strain background, and in which locus an epitope tag was added')
	parser.add_argument('-i','--input-dir', metavar='input_dir', required=True, help='the directory where all the EpitopeID output files were saved (*-ID.tab)')
	
	args = parser.parse_args()
	return(args)

#ENCFF000ZHH	-	-	/libraries/ENCLB119JWQ/	02b95214b189b0ece1b8ad292d664e16	2010-10-20T00:00:00.000000+00:00	3115882	-	revoked	/biosamples/ENCBS622ENC/	/genetic-modifications/ENCGM980FIG/biosample-types/cell_line_EFO_0001200/	MCF 10A	/targets/SRC-human/	-	insertion	tagging	ER	/labs/kevin-struhl/	To generate v-Src:ER, DNA sequences encoding the Schmidt Ruppin A form of v-Src (a gift from Josh Kaplan) were fused to hbER and cloned into the pWZLblast3, pBabepuro3, or pLNCX retroviral expression vectors. In the case of pWZLblast3 (a gift from Jay Morgenstern), expression of a bicistronic mRNA from the murine leukemia virus long terminal repeat containing the encephalomyocarditis virus internal ribosome entry site permits the expression of both v-Src:ER and the blasticidin drug resistance gene - from PMID:9891045.
#ENCFF747KWI	2	/files/ENCFF169UAL/	/libraries/ENCLB778MLD/	01e1770124ed4d81b69c5b0a57a9f6c8	2016-04-29T01:33:57.351935+00:00	53859571	101	released	/biosamples/ENCBS237PTC/	/genetic-modifications/ENCGM317SDT/	/biosample-types/cell_line_EFO_0001182/	HEK293	/targets/HIC1-human/	site-specific recombination	insertion	tagging	eGFP	/labs/michael-snyder/	This construct contains an eGFP N-terminal tag to HIC1 under control of a CMV promoter, introduced to expression vector pcDNA5/FRT as part of Thermo Fisher's Flp-In system. It is stably integrated at the introduced FLP recombinase target (FRT) site in the host genome.
#ENCFF093AOC	1	/files/ENCFF146YVN/	/libraries/ENCLB013WRU/	121eaff8748c0db6c2e7a7d7433c6bf5	2016-04-28T09:31:32.704210+00:00	30343274	101	released	/biosamples/ENCBS045SQX/	/genetic-modifications/ENCGM426QTG/	/biosample-types/cell_line_EFO_0001182/	HEK293	/targets/ZSCAN16-human/	site-specific recombination	insertion	tagging	eGFP	/labs/michael-snyder/	This construct contains an eGFP N-terminal tag to ZSCAN16 under control of a CMV promoter, introduced to expression vector pcDNA5/FRT as part of Thermo Fisher's Flp-In system. It is stably integrated at the introduced FLP recombinase target (FRT) site in the host genome.
if __name__ == "__main__":
	'''Collect metadata and EpitopeID results to get detection stats on the eGFP ENCODE data'''
	args = getParams()
	dbtag2mtag = {'LAP-tag':'eGFP'}

	# Initialize tracking variables
	total_se_pe = [0,0]
	se_epitope_found_count = 0
	pe_epitope_found_count = 0
	pe_correct_target_count = 0
	# Parse metadata
	reader = open(args.metadata, 'r', encoding='utf-8')
	for mline in reader:
		mtokens = mline.strip().split('\t')
		# Skip paired-end Read 2 entries (PE captured on Read 1)
		if(mtokens[1]=='2'):
			continue
		
		# Pull relevant info from metadata tokens
		encff = mtokens[0]
		target = mtokens[14].split('/')[2].split("-")[0].upper()
		tags = mtokens[18].split('|')
		
		# Only report eGFP
		if("eGFP" not in tags):
			continue
		
		# Check file exists
		id_file = join(args.input_dir,"%s_R1-ID.tab" % encff)
		if(not isfile(id_file)):
			sys.stderr.write("%s: no results generated.\n" % (id_file))
			continue
		
		# Initialize id file variables to save
		id_read_counts = []
		target_found = 0
		matched_target_strings = []
		
		# Parse id file for target info & epitope read count
		id_reader = open(id_file,'r')
		line = id_reader.readline()
		while( line!="" ):
			if( line.find("EpitopeID\tEpitopeCount")==0 ):
				line = id_reader.readline()
				tokens = line.strip().split('\t')
				mtag2count = {}
				# Add count values parsed from ID results file
				while( len(tokens)==2 ):
					mtag = dbtag2mtag.get(tokens[0],tokens[0])
					count = int(tokens[1])
					mtag2count.update({mtag:count})
					# Iterate to next line
					line = id_reader.readline()
					tokens = line.strip().split('\t')
				# Set-up default counts of all expected tags to zero
				for t in tags:
					mtag2count.setdefault(t,0)
				# Turn tagcountinfo into an ordered list of strings
				for t in sorted(mtag2count,key=mtag2count.get):
					id_read_counts.append("%s:%i" % (t,mtag2count[t]))
				# dictionaries preserve insertion order in Python 3.7+
				# recreate dictionary on order of value size
				#{k: v for k, v in sorted(mtag2count.items(), key=lambda item: item[1])}
			
			if( line.find("GeneID\tEpitopeID\tEpitopeLocation")==0 ):
				line = id_reader.readline()
				tokens = line.strip().split('\t')
				rank = 1
				while( len(tokens)==5 ):
					if( tokens[0].upper().find(target)==0 ):
						if( target_found==0 ):
							target_found = rank
						matched_target_strings.append(tokens[0].split("|")[0])
					elif( tokens[0]!="GeneID" ):
						sys.stderr.write('%s:missed_target=%s, target=%s\n' % (encff,tokens[0], target) )
					line = id_reader.readline()
					tokens = line.strip().split('\t')
					rank += 1
			line = id_reader.readline()
		id_reader.close()
		# Output: EID\t[SE|PE]\tEpitopeReadCount\t[Target\tTargetsFound
		sys.stderr.write("|".join(matched_target_strings))
		sys.stdout.write( "\t".join([encff,mtokens[1],target,str(target_found),"|".join(id_read_counts),"\t".join(mtokens[3:])]) + "\n")
#		sys.stdout.write('%s\t%s\t%i\t%s\t%i\n' % ( eid,mtokens[1],id_read_count,target,target_found  ) )
		# Update summary counts
#		if( mtokens[1]=='-' ):
#			total_se_pe[0] += 1
#			if( id_read_count>0 ):
#				se_epitope_found_count += 1
#		if( mtokens[1]=='1' ):
#			total_se_pe[1] += 1
#			if( id_read_count>0 ):
#				pe_epitope_found_count += 1
#			if( target_found==1 ):
#				pe_correct_target_count += 1
		id_reader.close()
	reader.close()
	# Output Summary
#	sys.stdout.write('\nOutput Summary:\n')
#	sys.stdout.write('%i of %i SE datasets found epitope\n' % (se_epitope_found_count,total_se_pe[0]))
#	sys.stdout.write('%i of %i PE datasets found epitope\n' % (pe_epitope_found_count,total_se_pe[1]))
#	sys.stdout.write('%i of %i PE datasets identified correct locus\n' % (pe_correct_target_count,total_se_pe[1]))
