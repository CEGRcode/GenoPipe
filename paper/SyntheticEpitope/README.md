# Simulate Paired-End datasets for EpitopeID and evaluate performance


## Generate Synthetic Genomes
Create synthetic genomes to simulate from by generating a random 500bp "epitope" sequence
Insert it into a genome adjacent to the coding sequence of a gene
-sacCer3 genome insert Rap1-Nterm
-sacCer3 genome insert Reb1-Cterm
-hg19 genome insert CTCF-Nterm
-hg19 genome insert POL2HR-Nterm

## Sequencing-Depth Test:
For each synthetic genome, simulate 1000 datasets for each of the various number of reads
-sacCer3 genomes sequence at 10K, 100K, 1M, and 10M reads
-hg19 genomes sequence at 100K, 1M, 10M, and 50M reads
Run them through EpitopeID and evaluate how often EpitopeID correctly identifies the
location of the inserted sequence. Time EpitopeID to determine how its execution speed
performance.

## Mixed/Contamination Simulations
For each organism, subsample the datasets generated above and mix them in various ratios.
-sacCer3 mix reads from Rap1-Nterm & Reb1-Cterm at a sequencing depth of
-hg19 mix reads from CTCF-Nterm & POL2HR-Nterm at a sequencing depth of
-each pair of sets should be mixed in the following ratios: (90-10%, 80-20%, .., 10-90%)
Run the new "contaminated" datasets through EpitopeID and evaluate how often EpitopeID
correctly identifies the location of the inserted sequence of each population.


## Job execution order
```
# Create synthetic genomes
bash job/generate_synthetic_genomes.sh

# Set-up job scripts for depth simulations
bash job/build_jobs.sh   # will create a bunch of PBS scripts in the job/ directory. Based on depth_template.pbs and epitopeid_template.pbs

# Run depth simulations to create FASTQ input files (yeast & human)
qsub run_depth_1_Reb1-Cterm_R500.pbs
qsub run_depth_2_Rap1-Nterm_R500.pbs
qsub run_depth_X_....
...

# Run EpitopeID on depth simulations to create the reports (yeast & human)
qsub run_EpitopeID_1_Reb1-Cterm_R500.pbs
qsub run_EpitopeID_2_Rap1-Nterm_R500.pbs
qsub run_EpitopeID_X_...
...

# Compile results from depth simulations
bash job/compile_results.sh

# Run mixture simulations to create FASTQ input files
qsub job/run_mix_yeast.pbs
qsub job/run_mix_human.pbs

# Run EpitopeID on mixture simulations to create the reports
qsub job/run_EpitopeID_on_mix_yeast.pbs
qsub job/run_EpitopeID_on_mix_human.pbs

# Compile results from mixture simulations
bash job/compile_mix_results.sh
```
