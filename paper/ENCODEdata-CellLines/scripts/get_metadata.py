import json, requests, sys

def search_url(url):
	# Force return from the server in JSON format
	headers = {'accept': 'application/json'}

	# GET the search result
	response = requests.get(url, headers=headers)

	# Extract the JSON response as a python dictionary
	search_results = response.json()
	return(search_results)

if __name__=="__main__":
	'''Retrieve metadata from ENCODE for testing StrainID'''

	uniq_cell_lines = []
	cell_lines2keep = ['A549','HCT116','HeLa','HeLa-S3','HepG2','K562','LNCAP','MCF-7','SK-N-SH']

	# Get Map for Biosample Types to Cell Line name
	sys.stderr.write("URL retrieve Biosample Types to create a map to string...")
	bst_url = "https://www.encodeproject.org/search/?type=BiosampleType&classification=cell+line&frame=object&format=json&limit=all"
	bst_search_results = search_url(bst_url)

	bst2name = {}
	sys.stderr.write("Searching %i Biosample Types found.\n" % len(bst_search_results["@graph"]))
	for bs_type in bst_search_results["@graph"]:
		bst_id = bs_type["@id"]
		term_name = bs_type["term_name"]
		if(term_name not in uniq_cell_lines):
			uniq_cell_lines.append(term_name)
		if(term_name in cell_lines2keep):
			bst2name.update({bst_id:term_name})
	bst2keep = bst2name.keys()

	# This searches the ENCODE database for the File codes with the file_format="bam" and assembly="hg19" criteria
	sys.stderr.write("URL retrieve all BAM files based on hg19 assembly...\n")
	ff_url = "https://www.encodeproject.org/search/?type=File&file_format=bam&assembly=hg19&frame=object&format=json&limit=all"
	ff_search_results = search_url(ff_url)

	sys.stderr.write("Searching %i BAM files found.\n" % len(ff_search_results["@graph"]))
	for bam in ff_search_results["@graph"]:
		#print(json.dumps(fastq, indent=4))
		accession = bam["accession"]

		# Collect File object info and store to a string (skip HiC assays)
		href = bam.get("href","-")
		md5sum = bam.get("md5sum","-")
		date_created = bam.get("date_created","-")
		datset = bam.get("dataset","-")
		biosample_ontology = bam.get("biosample_ontology","-")
		assay_term_name = bam.get("assay_term_name","-")
		output_type = bam.get("output_type","-")
		# read_count = str(bam.get("read_count","-"))
		file_size = bam.get("file_size","-")
		mapped_run_type = bam.get("mapped_run_type","-")
		mapped_read_length = bam.get("mapped_read_length","-")
		status = bam.get("status","no_status")
		derived_from = "|".join(bam.get("derived_from",[]))


		if(isinstance(assay_term_name,list)):
			assay_term_name = "|".join(assay_term_name)

		#if(assay_term_name != "ChIP-seq"):
		#	continue

		if(isinstance(biosample_ontology,list)):
			cell_line = "|".join([bst2name.get(i,i) for i in biosample_ontology])
			biosample_ontology = "|".join(biosample_ontology)
		elif(biosample_ontology not in bst2keep):
			continue
		else:
			cell_line = bst2name[biosample_ontology]

		info_tokens = [ href, md5sum, cell_line, assay_term_name, mapped_run_type, mapped_read_length, derived_from, output_type,  datset, biosample_ontology, file_size, status, date_created ]
		ff_info = "\t".join([ str(i) for i in info_tokens ])

		# Write info and accessions to a metdata data entry line
		sys.stdout.write( "\t".join([accession,ff_info]) + "\n" )
