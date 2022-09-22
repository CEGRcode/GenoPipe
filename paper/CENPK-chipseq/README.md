# StrainID detection of background strain in CEN.PK ChIP-seq samples

Run StrainID on CENPK datasets to evaluate StrainID's ability to detect the variant-based strain background.

["Integration of multiple nutrient cues and regulation of lifespan by ribosomal transcription factor Ifh1" (Cai et al, 2013)](https://pubmed.ncbi.nlm.nih.gov/24035395/)
### GEO ID : [GSE39147](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39147)


## Download data using SRA accessions using `job/00_download_data.pbs`

## Align FASTQ files and process using `job/01_align_data.pbs`

## Run StrainID on BAM inputs using `job/02_run_StrainID.pbs` to determine if StrainID can successfully identify the strain background
- The default sacCer3 StrainID database is used
