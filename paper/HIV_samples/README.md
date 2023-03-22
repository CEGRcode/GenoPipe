# EpitopeID detection of HIV genome insertions

Run EpitopeID on HIV datasets to evaluate EpitopeID's ability to detect HIV genome insertions.

["Benzotriazoles Reactivate Latent HIV-1 through Inactivation of STAT5 SUMOylation" (Bosque et al, 2017)](https://pubmed.ncbi.nlm.nih.gov/28147284/)

### GEO ID : [GSE84199](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84199)

## EpitopeID database with HIV genome as a tag and hg19 genome as genomic sequence is setup with `../setup.sh`
[Download the HIV genome (AF324493)...](https://www.ncbi.nlm.nih.gov/nuccore/AF324493)

## Download data using SRA accessions using `job/00_download_data.pbs`

## Run EpitopeID on FASTQ inputs using `job/01_run_EpitopeID.pbs` to determine if EpitopeID can localize HIV genome insertions
