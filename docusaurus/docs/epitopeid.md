---
id: epitopeid
title: EpitopeID
sidebar_label: EpitopeID
---

<!-- ![epitopeid-icon](/genopipe-img/epitopeid-icon.png)-->

__EpitopeID identifies and determines the genomic location of epitopes or other inserted sequences relative to genomic loci.__

![Figure1A](/genopipe-img/figure1a.png)

## Usage
```bash
bash identify-Epitope.sh -i /path/to/FASTQ -o /path/to/output -d /path/to/genome/database -t <Threads - Default 1>
eg: bash identify-Epitope.sh -i /input -o /output -d /sacCer3_EpiID -t 2
```


Specific Dependencies


## Input (`-i`)

EpitopeID takes [gzipped][gzip-man] FASTQ files from single-end(SE) or paired-end(PE) datasets to run. These should all be added to the same directory and the path to this directory will be used as the input for `identify-Epitope.sh`. If your FASTQ files are not already gzipped and you have `gzip` installed, you can simply gzip them in your terminal:

```bash
gzip XXXXX_R1.fastq
```

It is expected that at least one file for each sample(even SE) follows the naming convention of `_R1*.fastq.gz` and a second file if the data is paired-end following `_R2*.fastq.gz`. This is based on the Illumina naming standard.

:::caution
Make sure none of the sample names used in the filenames have an occurrence of `_R1` to avoid errors from EpitopeID about being unable to find files or attempting to use a `_R2` file as a `_R1` file.
:::



## Reference Files (`-d`)

:::note
The provided database files are missing the genomic reference file. You will need to follow the directions below to download `genome.fa` before running EpitopeID if you are planning to use the provided references.
:::

For EpitopeID, this is the "database" or directory with four types of reference files used by `identify-Epitope.sh`. You will notice that EpitopeID provides reference files for both yeast (`sacCer3_EpiID`) and human (`hg19_EpiID`) so you can quickly get started without building up the database from scratch. However, you are free to customize and build your own set of files (e.g. add different epitope tags to check for, use a different genome build).

Below is a list of the the hardcoded filenames that EpitopeID looks for during execution and some information on the provided yeast and human defaults.

* The `FASTA_tag/ALL_TAG.fa` is the FASTA formatted collection of all epitope sequences to search for. The yeast tag database includes the [AID, Extended-Tap,  FRB, HA_v1, HA_v3, MNase, ProtA, CBP, FLAG-3x, GFP, HA_v2, HaloTag, and Myc(3x)][tag-ref] sequences. The human tag database only includes the [LAP][lap-ref] tag but it is easy to customize the list to include other epitopes for EpitopeID to look for.
* The `FASTA_genome/genome.fa` is the reference genome used for the genomic alignments that the other annotations are base on. Even if you use the provided databases from Github, the genomic reference file still needs to be downloaded and moved to `FASTA_genome/genome.fa`. (Genome was not include for data storage reasons)
* The `annotation/genome_annotation.gff.gz` file defines the bin coordinates to use when localizing the epitope insertion in PE datasets. The yeast default uses SGD gene annotation coordinates to defines one bin for the length of each gene, 250bp bins flanking each set of gene coordinates, and 250bp bins breaking up the remaining intergenic regions. The human default similarly bins out the genome using 1000bp windows on the NCBI Refseq annotations.
* The `blacklist_filter/blacklist.bed`


Whether you use the provided reference files or create your own, the database should use the following directory structure both to ensure that EpitopeID can find the correct reference files and for organization, clarity, and consistency.
```
/name/of/epiDB
|--FASTA_tag
|  |--Tag_DB
|  |  |--tagname1.fa
|  |  |--tagname2.fa
|  |  |--tagname3.fa
|  |--ALL_TAG.fa
|  |--ALL_TAG.fa.amb   # EpitopeID will automatically generate
|  |--ALL_TAG.fa.ann   # these BWA index files if any are missing
|  |--ALL_TAG.fa.bwt   #
|  |--ALL_TAG.fa.fai   #
|  |--ALL_TAG.fa.pac   #
|  |--ALL_TAG.fa.sa    #
|--FASTA_genome
|  |--genome.fa
|--annotation
|  |--genome_annotation.gff.gz
|--blacklist_filter
|  |--blacklist.bed
```

Below is more information on how to use the utility scripts to download and customize your reference files.


### Downloading the sacCer3 genome

Use the utility scripts by following the commands below to download the sacCer3 genome and format the chromosome names. The reference files are based on the arabic numeral chromosome naming system (i.e. "chr1", "chr2",..."chr16","chrM","2-micron").

```bash
cd /path/to/GenoPipe/EpitopeID/utility_scripts/genome_data
bash download_sacCer3_Genome.sh
mv genome.fa /path/to/sacCer3_EpiID/FASTA_genome/
```

### Downloading the hg19 genome

Use the utility scripts by following the commands below to download the hg19 genome and move it to the appropriate directory.

```bash
cd /path/to/GenoPipe/EpitopeID/utility_scripts/genome_data
bash download_hg19_Genome.sh
mv genome.fa /path/to/hg19_EpiID/FASTA_genome/genome.fa
```


### Customizing epitopes

If you choose to add or remove epitope tags to your database, you must add/remove the files with the sequences in FASTA format to/from `FASTA_tag/TagDB` and recreate the `FASTA_tag/ALL_TAG.fa` file so it includes the tag sequences you want. The following commands show you how to remake the file using the provided scripts.

```bash
# Copy or remove the FASTA files with your tag sequences into the Tag_DB directory
cp tag1.fa /path/to/Tag_DB
cd /path/to/Tag_DB
bash /path/to/GenoPipe/EpitopeID/utility_scripts/update_TagDB.sh
mv ALL_TAG.fa* /path/to/hg19_EpiID/FASTA_tag/
```



### Customizing annotation


#### Make `genome_annotation.gff.gz` with a different bin size


:::warning
write up command series for yeast
```bash
# bedtools intersect -wb -abam $OUTPUT/$SAMPLE/orf_filter.bam -b $DATABASE/annotation/genome_annotation.gff.gz -bed > $OUTPUT/$SAMPLE/align-pe.out
cd /path/to/GenoPipe/EpitopeID/utility_scripts/annotation_data
bash

```
:::



:::warning
write up command series for human
:::

#### Make `genome_annotation.gff.gz` with a different set of annotations


:::warning
write up command series for unicorn
:::

```bash
# bedtools intersect -wb -abam $OUTPUT/$SAMPLE/orf_filter.bam -b $DATABASE/annotation/genome_annotation.gff.gz -bed > $OUTPUT/$SAMPLE/align-pe.out
cd /path/to/GenoPipe/EpitopeID/utility_scripts/annotation_data
bash

```



### Customizing filter

```bash
# bedtools intersect -v -abam $OUTPUT/$SAMPLE/orf.bam -b $DATABASE/blacklist_filter/blacklist.bed > $OUTPUT/$SAMPLE/orf_filter.bam
mv /path/to/blacklist_filter.bed /path/to/EpiID-DB/blacklist_filter/blacklist.bed
```



<!-- ### Custom genome from scratch

```bash
EPITOPEID=/path/to/GenoPipe/EpitopeID
# outline new directory structure
mkdir unicorn_EpiID
cd unicorn_EpiID
mkdir FASTA_genome
mkdir -p FASTA_tag/TagDB
mkdir annotation
mkdir blacklist_filter
# copy over genome
cp /path/to/unicorn/genome.fa FASTA_genome/genome.fa
# copy over tag sequences
cp /path/to/tag/sequences/*.fa FASTA_tag/TagDB/
cd $EPITOPEID/utility_scripts

``` -->

## Example: Set-up EpitopeID and run on yeast example
```bash
git clone www.github/CEGRcode/GenoPipe
cd GenoPipe/EpitopeID/

# remake the epitope Tags file and index it
cd sacCer3_EpiID/FASTA_tag/Tag_DB
bash ../../../utility_scripts/update_TagDB.sh
mv ALL_TAG.fa* ../
# download the genome
cd ../../../utility_scripts/genome_data/
bash download_sacCer3_Genome.sh
mv genome.fa ../../sacCer3_EpiID/FASTA_genome/
# download annotations
# cd ../annotation_data/
# bash generate_sacCer3_GenomeAnnotation.sh -g ../../sacCer3_EpiID/FASTA_genome/genome.fa -b 250
# mv genome_annotation.gff.gz ../../sacCer3_EpiID/annotation/
cd ../../
# run EpitopeID on FASTQs in the `sample` directory below
mkdir ../output
bash identify_Epitope.sh -i ../samples/ -o ../output/ -d sacCer3_EpiID -t 4
```






## Threading (`-t`)

This optional input is used to specify the number of threads to used for the BWA alignment commands.



## Output Report (`-o`)

The output report is saved to the user-provided output directory in a file named based on the input FASTQ files (`/path/to/output/samplename_R1-ID.tab`). Below is a sample report based on the results from running EpitopeID on the ENCODE ENCFF415CJF sample.

```
EpitopeID	EpitopeCount
LAP-tag	435

GeneID	EpitopeID	EpitopeLocation	EpitopeCount	pVal
NR4A1|chr12:52416616-52453291	LAP-tag	C-term	9	3.580493355965414e-24
```

The first part of the report shows which epitopes in `Tag_DB` were identified in the sample (**EpitopeID column**) and how many reads mapped to this epitope (**EpitopeCount**) to help quantify the coverage of the epitopes which relates to the confidence of the call.

The second part of the report shows which epitopes localized to which regions/tiles of the genome significantly (sorted by pvalue if multiple hits). The columns specify the cooridinate interval (**GeneID**), which epitope maps to this locus (**EpitopeID**), if this occurs on the N or C-terminus (**EpitopeLocation**), the number of reads mapping to this tile (**EpitopeCount**), and the poisson-calculated associated p-value to indicate confidence of the site (**pVal**).



## FAQs

* Q: I added my own custom tag sequences to the `TagDB` directory but when I run EpitopeID, none of my samples are getting significant hits to the new tags.
  * There are a few things you should check before concluding that the epitope is not present in your sample:
  * Did you recreate the `ALL_TAG.fa` file? Open it up to make sure your sequences are there. If they aren't there, follow the commands in the "How to add your own epitope tag sequences" section above.
  * Did you recreate the BWA index files for `ALL_TAG.fa`? Remove at least one of the BWA index files (e.g. `ALL_TAG.fa.ann`) and run EpitopeID. EpitopeID automatically recreates the index if any files are missing but it does not check that the index files match the FASTA sequence. If you modify the `ALL_TAG.fa` without recreating the BWA index files, EpitopeID will run on the old set of tag sequences.
* Q: What does the N-term or C-term mean in the output?
  * For each gene, three bins are created for (1) the ORF interval, (2) a bin upstream of the start codon, (3) a bin downstream of the stop codon, each of which correspond to the gene, N-terminus, or C-terminus of the peptide chain gene product. If EpitopeID is mapping more strongly to the N or C-terminus, it is likely that the epitope is tagged to that side of the endogenous peptide chain.
* Q: I tried running EpitopeID on my samples using `sacCer3_EpiID`/`hg19_EpiID` and it is telling me that I am missing the `genomic.fa` file. Where can I find it?
  * Due to storage reasons, we do not include the genomic reference in the Github download for GenoPipe. However, we do provide scripts that will help you to download and format the reference genomes that we used. See above for directions on how to use them by looking at the section titled "Downloading the sacCer3/hg19 genome". Otherwise you can use your own genomic sequence. If you do, then remember to rename it to `genome.fa` and move it to the appropriate directory. Also make sure that the annotations you use are based on the same reference build.


[gzip-man]:https://www.gnu.org/software/gzip/manual/gzip.html
[AID-ref]:https://www.google.com
[Extended-Tap-ref]:https://www.google.com
[FRB-ref]:https://www.google.com
[HA_v1-ref]:https://www.google.com
[HA_v3-ref]:https://www.google.com
[MNase-ref]:https://www.google.com
[ProtA-ref]:https://www.google.com
[CBP-ref]:https://www.google.com
[FLAG-3x-ref]:https://www.google.com
[GFP-ref]:https://www.google.com
[HA_v2-ref]:https://www.google.com
[HaloTag-ref]:https://www.google.com
[Myc(3x)-ref]:https://www.google.com
[tag-ref]:https://www.google.com
[lap-ref]:https://www.google.com
