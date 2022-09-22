# Run EpitopeID on ENCODE data and evaluate EptiopeID's performance


# Reference files
ENCODE metadata was pulled on February 27, 2021 using the `scripts/get_metadata.py`
script that pulls all Genetic Modification accessions with category=insertion and 
purpose=tagging. These are used to pull Biosample accessions with the organism=
human constraint. The File accessions with type=FASTQ are finally pulled with the 
filter that they are associated with Libraries that come from one of these Biosample
accessions and their assay=[ChIPseq|DNAseq]. They are saved with all relevant 
metadata to the `210227_sample_metadata.txt` file.

Command used: python scripts/get_metadata.py > 210227_sample_metadata.txt

# Download ENCODE eGFP data and run through EpitopeID
Use the files.txt with ENCODE accessions to download the data to the data directory.
Then make sure the LAP-tag is in the hg19_EpiID and run the data through EpitopeID.

# Compare the results to the metadata information
Evaluate the accuracy of EpitopeID on real data.
