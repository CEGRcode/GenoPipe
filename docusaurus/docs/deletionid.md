---
id: deletionid
title: DeletionID
sidebar_label: DeletionID
---

<!-- ![deletionid-icon] -->

__DeletionID identifies significant depletions of aligned NGS tags in the genome relative to a background model. This tool is set-up to confirm gene knockouts.__

![Figure1B]

## Usage
```bash
bash identify-Deletion.sh -i /path/to/BAM -o /path/to/output -d /path/to/genome/database
eg: bash identify-Deletion.sh -i /input -o /output -d /sacCer3_Del
```


## Input(`-i`)
DeletionID takes a directory pathname for the input and will run DeletonID on all the BAM files (`*.bam`) it finds within that directory path when you execute `identify-Deletion.sh`.

:::caution
Using __RNAseq__ data are not recommended!

Because DeletionID works by looking for significant depletion of reads in the genome, NGS assays with read distributions across the genome that dramatically deviate from uniform coverage will not work well. However assays like ChIP-exo for which we have validated this tool on are appropriate inputs. This is due to the large background of ChIP-based assays and the smoothing effect of binning the strategy that DeletionID implements. (Will?) Generally, genomic sequencing assays work well and transcriptomic assays are to be avoided for this module.
:::

:::caution
Identifying deletions of highly __repetitive regions__ of the genome using DeletionID is not recommended!

These highly repetitive regions of the genome have very low mappability scores and are often thrown out without even trying to calculate the depletion score.
:::

:::caution
If you are trying to run DeltionID on samples from __human__ or another organism with a very large genome, please read this!

DeletionID requires a higher amount of sequencing coverage of the genome to get reliable results. Because of both this and the very small size of human gene knockout mutations (esp relative to the size of the genome), everything might come back  with a significant you'll  likely just be measuring a bunch of noise.
:::




## Reference Files (`-d`)

For DeletionID, this is the "database" or directory with all the reference files used by `identify-Deletion.sh`. You will notice that DeletionID provides reference files for both yeast (`sacCer3_Del`) and human (`hg19_Del`) so you can quickly get started without building up the database from scratch. However, you are free to customize the database by adding new mappability scores or by looking at a different set of coordinate intervals.

Below is a list of the the hardcoded filenames that DeletionID looks for during execution and some information on the provided yeast and human defaults.

* The `coord.bed` is the file with the genomic coordinate intervals we are looking for depletions within. The provided `sacCer3_Del` uses the set of yeast gene coordinate intervals.
* The `mappability.tab` is the file with the calculated mappability scores for each coordinate interval within the `coord.bed` file. This will need to be regenerated if you use a different `coord.bed` file from the one provided.

Whether you use the provided reference files or create your own, the database should use the following directory structure both to ensure that DeletionID can find the correct reference files and for consistency.

```
/name/of/delDB
|--genomic_coord
   |--coord.bed
|--mappability
   |--mappability.tab
```

Below is more information on how to use the utility scripts to download and customize your reference files.

### How to add custom coordinate bed file

For custom coordinates, adding your BED-formatted genomic intervals to `yourDB/genomic_coord/coord.bed` file. Next you must generate new mappability scores for the new intervals. (See how to generate new mappabilities for a set of coordinates below)

### How to generate new mappabilities for a set of coordinates


(if you want to change the parameters for calculating mappability)

## Output Report (`-o`)



## FAQs

* Q: alsdjflasdkfj
* A: sequencing coverage is too low

* Q:
* A: DeletionID does not work well with RNAseq datasets.

* Q: I have data from a human/mouse/unicorn deletion
* A: DeletionID does not work well with samples from human or other organisms with large genomes. It has been tested and works fine in yeast for data 2M-10M PE reads but has not been tested in organisms with larger genomes. Please keep this in mind

* Q:
* A: does not work for high repeat regions---or maybe it does if you want to expand your deletion  to include the whole repeat region and lower the threshold of the Python script reporting to spit out all numbers and you can compare the score of the region to the scores of a control sample. This may require studying several replicates, perhaps even across several studies to determine if the numbers can be used to measure dramatic expansions...



[deletionid-icon]:../static/genopipe-img/deletionid-icon.png

[Figure1B]:../static/genopipe-img/figure1b.png
