# StrainID detection of background strain in BY4742 ChIP-seq samples

Run StrainID on BY4742 datasets to evaluate StrainID's ability to detect the variant-based strain background

[ "Molecular mechanisms that distinguish TFIID housekeeping from regulatable SAGA promoters" (de Jonge et al, 2017)](https://pubmed.ncbi.nlm.nih.gov/27979920/)
### GEO ID : [GSE81787](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81787)


## Download data using SRA accessions using `job/00_download_data.pbs`

## Align FASTQ files and process using `job/01_align_data.pbs`

## Run StrainID on BAM inputs using `job/02_run_StrainID.pbs` to determine if StrainID can successfully identify the strain background
- The default sacCer3 StrainID database is used
- specifically look into performance in distinguishing two closely related strains: BY4741 and BY4742
