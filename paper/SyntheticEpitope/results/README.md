# Simulation files and EpitopeID results go here

## sacCer3
This directory contains the raw EpitopeID reports and runtimes for each of the 1000 simulations across each type of yeast simulation.

## hg19
This directory contains the raw EpitopeID reports and runtimes for each of the 1000 simulations across each type of human simulation.

## Summary reports
The EpitopeID results are parsed out into a summary report using the `scripts/analyze_eid_results.py` where...
- **First Column:** is the filepath for each simulation's EpitopeID report (filepath encodes experiment parameters)
- **Second Column:** includes the read count of reads that map to the expected epitope sequence
- **Third Column:** includes the number of bins relating to the expected target region (>0 means successfully localized)


The following files contain the results from variable sequencing depth and epitope tag length experiments in yeast and human:
```
SummaryReport_sacCer3.txt
SummaryReport_hg19.txt
```
![sacCer3-id-tally](results/ID-tally_sacCer3.png)
![hg19-id-tally](results/ID-tally_hg19.png)

The following files contain the results from the read mixture contamination titration experiments in yeast and human:
```
MixSummaryReport_sacCer3.txt
MixSummaryReport_hg19.txt
```
![mix-sacCer3-id-tally](results/ID-Mix-tally_sacCer3_1M.png)
![mix-hg19-id-tally](results/ID-Mix-tally_hg19_50M.png)

## Runtime Summary reports
The following files contain the results from runtime benchmarking the read mixture contamination titration experiments in yeast and human:
```
RuntimeSummaryReport_sacCer3.txt
RuntimeSummaryReport_hg19.txt
```
![sacCer3-runtimes](results/Runtimes_sacCer3.png)
![hg19-runtimes](results/Runtimes_hg19.png)
