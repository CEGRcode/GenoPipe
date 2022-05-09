---
id: deletionid
title: DeletionID
sidebar_label: DeletionID
---

<!-- [deletionid-icon]:../static/genopipe-img/deletionid-icon.png -->

__DeletionID identifies significant depletions of aligned NGS tags in the genome relative to a background model. This tool is set-up to confirm gene knockouts.__

![Figure1B](/genopipe-img/figure1b.png)

## Usage
```bash
bash identify-Deletion.sh -i /path/to/BAM -o /path/to/output -d /path/to/genome/database
eg: bash identify-Deletion.sh -i /input -o /output -d /sacCer3_Del
```


## Input(`-i`)
DeletionID takes a directory pathname for the input and will run DeletionID on all the BAM files (`*.bam`) it finds within that directory path when you execute `identify-Deletion.sh`.

:::caution
Using __RNAseq__ data is not recommended!
:::

Because DeletionID works by looking for significant depletion of reads in the genome, NGS assays with read distributions across the genome that dramatically deviate from uniform coverage will not work well. However chromatin-based assays like whole genome sequencing (WGS) and ChIP-exo for which we have validated this tool on are appropriate inputs. This is due to the large background of ChIP-based assays and the smoothing effect of the binning strategy that DeletionID implements. Generally, genomic sequencing assays work well and transcriptomic assays are to be avoided for this module.

:::caution
If you are trying to run DeletionID on samples from __human__ or another organism with a larger genome than yeast, please read this!
:::

DeletionID has been tested and works fine in yeast for data >3M PE reads but has not been thoroughly tested in organisms with larger genomes. Larger genomes will require a higher amount of sequencing to cover the genome and establish a consistent background model for checking relative depletion. In addition, the small size of human gene knockout mutations relative to the size of the genome may result in poor identification of true deletions. Many intervals may be reported as significant that are not truly knocked out. You may also run into mappability problems due to the larger genomes tending to have more repetitive regions. Regions with low mappability may need to be blacklisted in the analysis. Please keep this in mind when running DeletionID and try to use high coverage datasets to confirm knockouts.


## Reference Files (`-d`)

For DeletionID, this is the "database" or directory with all the reference files used by `identify-Deletion.sh`. You will notice that DeletionID provides reference files for yeast (`sacCer3_Del`) so you can quickly get started without building up the database from scratch. However, you are free to customize the database by adding a different set of coordinate new mappability scores or by looking at a different set of coordinate intervals.

Below is a list of the files that DeletionID looks for during execution and some information on the provided yeast and human defaults.

* The `coord.bed` is the file with the genomic coordinate intervals we are looking for depletions within. The provided `sacCer3_Del` uses the set of yeast gene coordinate intervals.
* The `mappability.tab` is the file with the calculated mappability scores for each coordinate interval within the `coord.bed` file. This will need to be regenerated if you use a different `coord.bed` file from the one provided.

Whether you use the provided reference files or create your own, the database should use the following directory structure both to ensure that DeletionID can find the correct reference files.

```
/name/of/delDB
|--genomic_coord
   |--coord.bed
|--mappability
   |--mappability.tab
```

Below is more information on how to use the utility scripts to download and customize your reference files.

### How to add custom coordinate BED file

For custom coordinates, adding your BED-formatted genomic intervals to `/name/of/delDB/genomic_coord/coord.bed` file. Next you must generate new mappability scores for the new intervals. (See how to generate new mappabilities for a set of coordinates below)

:::caution
Identifying deletions of highly __repetitive regions__ of the genome using DeletionID is not recommended!
:::

These highly repetitive regions of the genome have very low mappability scores and are often thrown out even before calculating the depletion score. If you wish to proceed, please try to expand the interval to include the entire repeat region and lower the threshold of the Python script reporting to spit out all scores. You may want to compare the score of the region to the scores of a control sample. The results may be strengthened using replicates, perhaps even across several studies to determine if the numbers can be used to measure dramatic expansions...


### How to generate a new mappability reference file

For each set of coordinates, a companion mappability file must be generated for DeletionID. The following code shows how to generate this reference file using the coordinate interval file and the genome FASTA that the coordinates are based on.

```
cd /path/to/GenoPipe/DeletionID/utility_scripts
bash generate_Mappability_File.sh -f /path/to/genome.fa -c /path/to/delDB/genomic_coord/coord.bed -t <Threads>
mv mappability.tab /path/to/delDB/mappability/mappability.tab
```


## Output Report (`-o`)

The output report is saved to the user-provided output directory in a file named based on the input BAM files (`/path/to/output/XXXXX_deletion.tab` from some input `XXXXX.bam`). Below is a sample report based on the results from running DeletionID on a whole genome sequencing dataset of the APE3 deleted knockout from the Yeast Knockout Collection (ERS838258 sample downloaded from ENA and generated by Puddu et al, 2019).

```
LEU2 No Data Detected
URA3 No Data Detected
APE3 -3.3992739779627357
```

The first column of the report lists out the intervals (ORFs) with significantly depleted reads, sorted by the DeletionID score values in the second column which sort the more negative scores (more read depleted intervals) first. This sample is from an APE3 knockout in a LEU2 and URA3 knockout background so we expect all three gene intervals to show up in the report. LEU2 and URA3 show **No Data Detected** which indicates no reads mapped to this interval. APE3 shows a very strong depletion with a large negative score (-3.3992...) confirming the knockout as likely.


## FAQs

* Q: I have run DeletionID but many genes are being returned in the output report. Are all of these genes in the list depleted from my sample?
  * First check your sequencing coverage. When it is low, DeletionID can misidentify or identify many more knockout regions than are actually present.
  * Is your assay chromatin based? DeletionID does not work well with RNAseq datasets due to the heavily skewed distribution of reads across the genome. Since scores are normalized by the median coverage under the assumption that the median score sufficiently represents the typical background coverage of a random interval, RNAseq style coverage would falsely identify a significant number of knockout sites for the lowly-expressed genes.

* Q: I have data from a human/mouse deletion that I wish to identify. Do you have a reference database for genes from the mm10 or hg19 genome assemblies?
  * DeletionID does not work well with samples from human or other organisms with large genomes. Please see the caution note above for more information.
