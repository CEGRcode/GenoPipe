# Simulate Paired-End datasets for StrainID and evaluate performance

## Generate Synthetic Genomes
Create synthetic genomes to simulate from by using the VCFs from StrainID to replace ref
nucleotides with the variant nucleotides.
-sacCer3 genome insert CEN.PK2-1Ca alleles
-sacCer3 genome insert RM11-1A allelesa
-hg19 genome insert K562 alleles
-hg19 genome insert HCT116 alleles

## Sequencing-Depth Test:
For each synthetic genome, simulate 1000 datasets for each of the various number of reads
-sacCer3 genomes sequence at 500K, 1M, 2M, 3M, 4M, and 5M reads
-hg19 genomes sequence at 1M, 2M, 5M, 10M, and 20M reads
Run them through StrainID and evaluate how often StrainID correctly identifies the
correct strain. Time StrainID to determine how its execution speed performance.

## Mixed/Contamination Simulations
For each organism, subsample the datasets generated above and mix them in various ratios.
-sacCer3 mix reads from CEN.PK2-1Ca & RM11-1A at a sequencing depth of
-hg19 mix reads from K562 & HCT116 at a sequencing depth of
-each pair of sets should be mixed in the following ratios: (90-10%, 80-20%, .., 10-90%)
Run the new "contaminated" datasets through StrainID and evaluate how often StrainID
correctly identifies the background of each population.