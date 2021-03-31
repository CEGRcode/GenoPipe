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
	'''Retrieve metadata from ENCODE for testing EpitopeID'''
	# Get Map for Biosample Types to Cell Line name
	sys.stderr.write("URL retrieve Biosample Types to create a map to string...")
	bst_url = "https://www.encodeproject.org/search/?type=BiosampleType&frame=object&format=json&limit=all"
	bst_search_results = search_url(bst_url)
	
	bst2info = {}
	sys.stderr.write("Searching %i Biosample Types found.\n" % len(bst_search_results["@graph"]))
	for bs_type in bst_search_results["@graph"]:
		#print(json.dumps(bs_type, indent=4))
		bst_id = bs_type["@id"]
		term_name = bs_type["term_name"]
		bst2info.update({bst_id:term_name})
	
	# This searches the ENCODE database for the GeneticModification codes with the "insertion" category and "tagging" purpose
	sys.stderr.write("URL retrieve Genetic Modifications with the \"insertion\" category and \"tagging\" purpose...\n")
	gm_url = "https://www.encodeproject.org/search/?type=GeneticModification&category=insertion&purpose=tagging&frame=object&format=json&limit=all"
	gm_search_results = search_url(gm_url)
	
	gm2info = {}
	bs2gm = {}
	
	sys.stderr.write("Searching %i genetic modifications found.\n" % len(gm_search_results["@graph"]))
	for modification in gm_search_results["@graph"]:
		#print(json.dumps(modification, indent=4))
		mid = modification["@id"]
		#"summary": "Homo sapiens K562 cell line genetically modified (insertion) using CRISPR targeting ZBTB17"
		
		#Pull fields from GM to save in gm2info
		target = modification["modified_site_by_target_id"]
		method = modification.get("method","-")
		category = modification["category"]
		purpose = modification["purpose"]
		description = modification.get("description","-")
		introduced_tags = modification.get("introduced_tags",[{"name":"-"}])
		tag_name_list = []
		for tag in introduced_tags:
			tag_name_list.append(tag["name"])
		tag_names = "|".join(tag_name_list)
		lab = modification.get("lab","-")
		# Check organism is human (will need to figure out why only name is being loaded in JSON later..)
		#organism = target["organism"]  #/organisms/human/
		#For now hacking through with the following line:
		#organism = target.split("-")[-1].split("/")[0]
		#if( organism != "human" ):
		#	sys.stderr.write("Skipping %s. Non-human organism %s.\n" % (accession, organism))
		#	continue
		
		gm_info = "\t".join([target,method,category,purpose,tag_names,lab,description])
		gm2info.update({mid:gm_info})
		
		#List out Biosamples associated with these GM
		biosamples = modification["biosamples_modified"]
		for bid in biosamples:
			bs2gm.update({bid:mid})
	
	bs_list = bs2gm.keys()
	
	# This searches the ENCODE database for the Biosample codes with the organisms="Homo sapiens" criteria
	sys.stderr.write("URL retrieve Biosamples with the \"Homo sapiens\" organism...\n")
	bs_url = "https://www.encodeproject.org/search/?type=Biosample&organism.scientific_name=Homo+sapiens&frame=object&format=json&limit=all"
	bs_search_results = search_url(bs_url)
	
	bs2info = {} 
	
	sys.stderr.write("Searching %i Biosamples found.\n" % len(bs_search_results["@graph"]))
	for biosample in bs_search_results["@graph"]:
		bid = biosample["@id"]
		ontology = biosample.get("biosample_ontology","-")
		term_name = bst2info.get(ontology,"-")
		# Save relevant info on Biosamples of interest
		if bid in bs_list:
			bs2info.update({bid:"\t".join([ontology,term_name])})
	
	# replace bs_list with filtered set and dereference to clear memory
	bs_list = bs2info.keys()
	
	# This searches the ENCODE database for the Biosample codes with the organisms="Homo sapiens" criteria
	sys.stderr.write("URL retrieve all Libraries...\n")
	lb_url = "https://www.encodeproject.org/search/?type=Library&frame=object&format=json&limit=all"
	lb_search_results = search_url(lb_url)
	
	lb2bs = {}
	
	sys.stderr.write("Searching %i Libraries found.\n" % len(lb_search_results["@graph"]))
	for library in lb_search_results["@graph"]:
		#print(json.dumps(library, indent=4))
		lid = library["@id"]
		accession = library["accession"]
		bid = library.get("biosample","no_biosample")
		
		# Save Library objects associated with our Biosamples of interest
		if( bid in bs_list ):
			lb2bs.update({lid:bid})
	
	lb_list = lb2bs.keys()
	
	# This searches the ENCODE database for the File codes with the file_format="fastq" criteria
	sys.stderr.write("URL retrieve all FASTQ files...\n")
	ff_url = "https://www.encodeproject.org/search/?type=File&file_format=fastq&frame=object&format=json&limit=all"
	ff_search_results = search_url(ff_url)
	
	sys.stderr.write("Searching %i FASTQ files found.\n" % len(ff_search_results["@graph"]))
	for fastq in ff_search_results["@graph"]:
		#print(json.dumps(fastq, indent=4))
		accession = fastq["accession"]
		
		# Skip Files not related to the Libraries of interest
		lid = fastq.get("library","-")
		if( lid not in lb_list ):
			#print(accession," ",lid)
			continue
		
		# Collect File object info and store to a string (skip HiC assays)
		paired_end = fastq.get("paired_end","-")
		paired_with = fastq.get("paired_with","-")
		md5sum = fastq.get("md5sum","-")
		assay_term_name = fastq.get("assay_term_name")
		date_created = fastq.get("date_created","-")
		read_count = str(fastq.get("read_count","-"))
		read_length = str(fastq.get("read_length","-"))
		status = fastq.get("status","no_status")
		date_created = fastq.get("date_created","-")
		
		if(assay_term_name=="HiC"):
			continue
		
		ff_info = "\t".join([paired_end,paired_with,lid,md5sum,assay_term_name,date_created,read_count,read_length,status])
		
		# Get accessions and info mappings for this file
		bid = lb2bs[lid]
		mid = bs2gm[bid]
		bs_info = bs2info[bid]
		gm_info = gm2info[mid]
		
		# Write info and accessions to a metdata data entry line
		sys.stdout.write( "\t".join([accession,ff_info,bid,mid,bs_info,gm_info]) + "\n" )
