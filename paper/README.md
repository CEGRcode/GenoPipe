# Lang et al, 202X

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

## setup.sh
runs the scripts to download and format the yeast and human genomes and other reference files for aligning the data
also indexes the genomes for BWA

## scripts
contains the general scripts for setting up the simulations like downloading and parsing the reference genomes
also contains the general scripts that several of the higher directory scripts call

## input
where setup.sh puts the reference genome and the aligner indexes
also where other reference FASTA files are stored (i.e. R500.fa, 3xFLAG.fa, and the HIV genome sequence)

## db
where the input database directories are built by setup.sh with variation as appropriate for GenoPipe module, species, and epitope set

## SyntheticEpitope
contains the scripts and houses the results of simulations testing EpitopeID

## SyntheticStrain
contains the scripts and houses the results of simulations testing StrainID

## SyntheticDeletion
contains the scripts and houses the results of simulations testing DeletionID

## ENCODEdata-eGFP
contains the scripts and information for downloading ENCODE eGFP data to test EpitopeID

## ENCODEdata-CellLines
contains the scripts and information for downloading ENCODE transcription factor ChIP-seq data to test StrainID

## HIV_samples
contains the scripts and information for downloading, processing, and running EpitopeID on the Bosque et al, 2017 dataset for localizing HIV genome insertions

## YKOC-wgs
contains the scripts and information for downloading, processing, and running DeletionID on the Puddu et al, 2019 dataset for identifying deletions

## BY4742-chipseq
contains the scripts and information for downloading, processing, and running StrainID on the BAM files

## CENPK-chipseq
contains the scripts and information for downloading, processing, and running StrainID on the BAM files
