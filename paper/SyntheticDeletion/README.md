# Simulate Paired-End datasets for DeletionID and evaluate performance

## Generate Synthetic Genomes
Create synthetic genomes to simulate from by deleting the coding sequence of a single
gene.
-sacCer3 genome delete Rap1
-sacCer3 genome delete Reb1
-hg19 genome delete CTCF
-hg19 genome delete POL2HR

## Sequencing-Depth Test:
For each synthetic genome, simulate 1000 datasets for each of the various number of reads
-sacCer3 genomes sequence at 10K, 100K, 1M, and 10M reads
-hg19 genomes sequence at 100K, 1M, 10M, and 50M reads
Run them through DeletionID and evaluate how often DeletionID correctly identifies the
location of the deleted sequence. Time DeletionID to determine how its execution speed
performance.

## Mixed/Contamination Simulations
For each organism, subsample the datasets generated above and mix them in various ratios.
-sacCer3 mix reads from Rap1-del & Reb1-del at a sequencing depth of
-hg19 mix reads from CTCF-del & POL2HR-del at a sequencing depth of
-each pair of sets should be mixed in the following ratios: (90-10%, 80-20%, .., 10-90%)
Run the new "contaminated" datasets through DeletionID and evaluate how often DeletionID
correctly identifies the location of the inserted sequence of each population.
