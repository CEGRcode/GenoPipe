# Lang et al (2023) *Nucleic Acids Res*

This contains all the scripts used to generate the figures validating and benchmarking performance for GenoPipe as a tool, specifically the code for downloading data, preprocessing, and simulations.

### This directory includes the simulation and preprocessing code used in the GenoPipe publication. If *you are a user* trying to just run GenoPipe the tool, navigate to either the EpitopID, StrainID, or DeletionID directories!

```
paper
|--setup.sh
|--scripts
|--input
|--db
|--SyntheticEpitope
|--SyntheticStrain
|--SyntheticDeletion
|--ENCODEdata-eGFP
|--ENCODEdata-CellLines
|--HIV_samples
|--YKOC-wgs
|--BY4742-chipseq
|--CENPK-chipseq
```

These scripts were built to run on a linux server with a PBS job scheduler set up and some of the dependencies installed using some environmental modules and a conda environment for remaining dependencies. You may need to modify these scripts to account for different server setup and configurations.

See the [GenoPipe documentation](https://pughlab.mbg.cornell.edu/GenoPipe-docs/) for a list of dependencies needed to run these publication-associated scripts. In addition to these dependencies, you will also need to install the following:

* [seqtk](https://github.com/lh3/seqtk).
* sra-toolkit (fastq-dump)
* wget

## setup.sh
Runs the scripts to download and format the yeast and human genomes and other reference files for aligning the data
also indexes the genomes for BWA. This must be run before all other scripts to reproduce the publication figures.

## scripts
Contains the general scripts for setting up the simulations like downloading and parsing the reference genomes
also contains the general scripts that several of the higher directory scripts call.

## input
Where `setup.sh` puts the reference genome and the aligner indexes. Also where other reference FASTA files are stored
(i.e. R500.fa, 3xFLAG.fa, and the HIV genome sequence).

## db
Where the input database directories are built by `setup.sh` with variation as appropriate for GenoPipe module, species, and epitope set.

## SyntheticEpitope
Contains the scripts and the results of simulations testing EpitopeID. This also includes the mixed contamination simulation tests.

## SyntheticStrain
Contains the scripts and the results of simulations testing StrainID.

## SyntheticDeletion
Contains the scripts and the results of simulations testing DeletionID.

## ENCODEdata-eGFP
Contains the scripts and information for downloading ENCODE eGFP data to test EpitopeID.

## ENCODEdata-CellLines
Contains the scripts and information for downloading ENCODE transcription factor ChIP-seq data to test StrainID.

## HIV_samples
Contains the scripts and information for downloading, processing, and running EpitopeID on the Bosque et al, 2017 dataset for localizing HIV genome insertions.

## YKOC-wgs
Contains the scripts and information for downloading, processing, and running DeletionID on the Puddu et al, 2019 dataset for identifying deletions.

## BY4742-chipseq
Contains the scripts and information for downloading, processing, and running StrainID on the BAM files.

## CENPK-chipseq
cContains the scripts and information for downloading, processing, and running StrainID on the BAM files.
