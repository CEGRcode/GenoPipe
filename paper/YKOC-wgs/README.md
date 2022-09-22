# Run DeletionID on WGS data of the YKOC

## Metadata with strain information, EBI accessions, and gene name maps.
The `201213_Puddu_2019_STable1_del2ers.txt` file was pulled from Supplementary
Table 1 from Puddu et al 2019. Maps the EBI sample accessions (ERSXXXXXXX) to the
expected knockout strain and a sample ID used by the paper (SDXXXX). This table was
downloaded December 13, 2020.

The `210328_Puddu_2019_STable6_verifiedKO.txt` file was pulled from Supplementary
Table 6 from Puddu et al 2019. This file contains the list of verified KO samples
keyed on their SD id and including associated scores with notes.

The `210403_PRJEB27160_accessions.txt` metadata file was pulled to map EBI sample
accessions (ERSXXXXXXX) to the run accessions (ERRXXXXXXX). The file also contains
other sequencing information including ftp locations for FASTQ files, MD5 checksums,
read counts, and read lengths for each sample. The following url was accessed on
April 3, 2021 to download this metadata.

https://www.ebi.ac.uk/ena/browser/view/PRJEB27160

Note: Four ERS sample accessions referenced in the paper are missing from the EBI
metadata download:

ERS1041599
ERS969339
ERS969700
ERS969701

## Download & Align YKOC data
Use the EBI accessions to download the data to the `results/FASTQ` directory. The
`job/00_download_data.pbs` PBS script uses the `wget` command but for the paper, these
files were downloaded using Globus. The FASTQ files are aligned with BWA-MEM by
running the `job/01_align_fastq.pbs` PBS script.

## Run YKOC through DeletionID
The BAM files are fed to DeletionID using the `job/02_indexed_runDID.pbs` PBS script.

## Compare the results to the metadata information
Evaluate the accuracy of DeletionID on real data using the `job/03_tally_results.pbs`
PBS script which calls the `script/analyze_ykoc_results.py` Python script.
