# Run EpitopeID on HIV datasets to evaluate EpitopeID's ability to detect HIV genome insertions

# "Benzotriazoles Reactivate Latent HIV-1 through Inactivation of STAT5 SUMOylation"
# (Bosque et al, 2017)

# GEO accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84199
# HIV genome: https://www.ncbi.nlm.nih.gov/nuccore/AF324493

# EpitopeID database with HIV genome as a tag and hg19 genome as genomic sequence is setup with `../setup.sh`
# Download data using SRA accessions using `job/00_download_data.pbs`
# Run EpitopeID on FASTQ inputs using `job/01_run_EpitopeID.pbs` to determine if EpitopeID can localize HIV genome insertions
