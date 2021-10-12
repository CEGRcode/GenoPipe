# Run StrainID on ENCODE data and evaluate StrainID's performance

# Reference files
ENCODE metadata was pulled on May 12, 2021 using the `scripts/get_metadata.py`
script that pulls all Biosample accessions with classification="cell_line" and
whose string matches one of the cell lines we have in our hg19_VCF database.
These are used to pull File accessions with type=BAM and assembly=hg19. We did
not filter by assay for this analysis. They are saved with all relevant metadata
to the `210512_sample_metadata.txt` file.

Command used: python scripts/get_metadata.py > 210512_sample_metadata.txt

# Download ENCODE CellLine data and run through StrainID
Use the 210512_sample_metadata.txt file with ENCODE accessions to download the
data to the data directory.
Then run the data through StrainID.

# Compare the results to the metadata information
Evaluate the accuracy of StrainID on real data.
