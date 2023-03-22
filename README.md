# GenoPipe

Expanded Documentation at https://pughlab.mbg.cornell.edu/GenoPipe-docs/

## Toolkit for characterizing the genotype of NGS datasets

There are 3 primary modules for genotype identification:

### EpitopeID

- Identify and determine the genomic location of epitopes relative to genomic loci

- Some epitope sequences are provided in the default tag database

| sacCer3(yeast) | hg19(human) |
| -------------- | ----------- |
| AID | LAP-tag |
| CBP |  |
| Extended-Tap |  |
| FLAG-3x |  |
| FRB |  |
| GFP |  |
| HA_v1 |  |
| HA_v2 |  |
| HA_v3 |  |
| HaloTag |  |
| MNase_v2 |  |
| Myc-3x |  |
| ProteinA |  |

### DeletionID

- Identify signficant depletion of aligned NGS tags in the genome relative to a background model. This module is useful for confirming gene knockouts.

- Default database includes reference files that direct the search for depleted reads within gene annotation intervals from the sacCer3(yeast) genome build.

### StrainID

- Compare a database of VCF files against an aligned BAM file to check for the presence of SNPs in order to determine a likely cell line/strain used in the experiment

- Default database includes reference VCF files for the following strains:

| sacCer3(yeast) | hg19(human) |
| -------------- | ----------- |
| BY4741         | A549        |
| BY4742         | HCT116      |
| CEN.PK2-1Ca    | HELA       |
| D273-10B       | HepG2        |
| FL100          | K562       |
| JK9-3d         | LnCap        |
| RM11-1A        | MCF7       |
| RedStar        | SKnSH        |
| SEY6210        |  |
| Sigma1278b-10560-6B |  |
| W303           |  |
| Y55            |  |


[Figure 1]


## Quickstart

This guide is for how to run each of the three GenoPipe modules on data from yeast(sacCer3) and human(hg19) samples. See the full documentation for how to modify and generate reference files for other genome builds.

### Dependencies

You will need the following software to run all GenoPipe modules:

[Samtools v1.5+](http://www.htslib.org/)

[Bedtools v2.27+](https://bedtools.readthedocs.io/en/latest/)

[BWA v0.7.15+](http://bio-bwa.sourceforge.net/bwa.shtml)

[Python v3.6.8+](https://www.python.org/)

- [scipy v1.5.4+](https://www.scipy.org/)

- [pysam v0.16.0.1+](https://pysam.readthedocs.io/en/latest/api.html)

[Perl](https://www.perl.org/)

[wget](https://www.gnu.org/software/wget/)

conda install:

```
conda create -n genopipe -c conda-forge -c bioconda python perl bwa bedtools samtools pysam scipy wget
```

### Download

To download GenoPipe, you can clone the repostitory. No builds needed.

```
git clone https://github.com/CEGRcode/GenoPipe
cd GenoPipe
```


###EpitopeID
genomes and epitope sequences

- yeast epitope tags

Saccer33xMyc



1. Check FASTQ filenames

EpitopeID takes gzipped FASTQ files as input. The file name should end with a `_R1` or `_R2` and use the extension `fastq.gz` (the standard naming convention of Illumina libraries).

Example:

The following would be valid file names for EpitopeID input files where “SampleA” is single-end data and SampleB is paired-ended.

```
SampleA_R1.fastq.gz
SampleB_R1.fastq.gz
SampleB_R2.fastq.gz
```

2. Set-up the database

The following instructions are for setting up the database of reference files used by EpitopeID using the provided genome builds and epitope tag sequences. To customize your database, see the full documentation.

For downloading yeast genome...

```
cd EpitopeID/utility_scripts/genome_data/
bash download_sacCer3_Genome.sh
mv genome.fa* ../../sacCer3_EpiID/FASTA_genome/
```

For downloading human genome...

```
cd EpitopeID/utility_scripts/genome_data/
bash download_hg19_Genome.sh
mv genome.fa* ../../hg19_EpiID/FASTA_genome/
```


3. Run EpitopeID

When providing path locations, it is important that you provide **absolute paths** (i.e. path should start with `/` or `~/`).

For yeast (sacCer3) samples...
```
cd GenoPipe/EpitopeID
bash identify-Epitope.sh -i /path/to/FASTQ -o /path/to/output -d /path/to/GenoPipe/EpitopeID/sacCer3_EpiID
```

For human (hg19) samples...
```
cd GenoPipe/EpitopeID
bash identify-Epitope.sh -i /path/to/FASTQ -o /path/to/output -d /path/to/GenoPipe/EpitopeID/hg19_EpiID
```


Joe Schmoe Example:

In the following example, GenoPipe, the directory including all the input yeast FASTQ files, and the new directory for storing EpitopeID reports are stored on the Desktop of Joe Schmoe. Filepaths would need to be changed according to a user's preferred directory structure.

```
# Download GenoPipe
cd /User/joeschmoe/Desktop/
git clone GenoPipe
# Download Genomic FASTA and move to appropriate directory
cd /User/joeschmoe/Desktop/GenoPipe/EpitopeID/utility_scripts/genome_data/
bash download_sacCer3_Genome.sh
mv genome.fa* ../../sacCer3_EpiID/FASTA_genome/
cd ../../
# Run EpitopeID
bash identify-Epitope.sh -i /User/joeschmoe/Desktop/myfastq -o /User/joeschmoe/Desktop/myreports_EID -d /User/joeschmoe/Desktop/GenoPipe/EpitopeID/sacCer3_EpiID
```




### DeletionID

1. Align FASTQ input files

DeletionID uses BAM files as its input. Make sure that the reads are aligned to sacCer3 if you are using the default interval database. Any aligner that outputs standard BAM format can be used to generate the BAM input. DeletionID was tested on [BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml).

2. Run DeletionID


For yeast (sacCer3) samples...

```
cd GenoPipe/DeletionID
bash identify-Deletion.sh -i /path/to/BAM -o /path/to/output -d /path/to/GenoPipe/DeletionID/sacCer3_Del
```

Joe Schmoe Example:

In the following example, GenoPipe, the directory including all the input yeast BAM files, and the new directory for storing DeletionID reports are stored on the Desktop of Joe Schmoe. Filepaths would need to be changed according to a user's preferred directory structure.

```
cd /User/joeschmoe/Desktop/GenoPipe/DeletionID
# Run DeletionID
bash identify-Deletion.sh -i /User/joeschmoe/Desktop/mybam -o /User/joeschmoe/Desktop/myreports_DID -d /User/joeschmoe/Desktop/GenoPipe/DeletionID/sacCer3_Del
```

### StrainID

1. Align FASTQ input files

StrainID uses BAM files as its input. Make sure that the reads are aligned to the appropriate sacCer3 or hg19 genome build if you are using the default interval database. Any aligner that outputs standard BAM format can be used to generate the BAM input. StrainID was tested on [BWA-MEM](http://bio-bwa.sourceforge.net/bwa.shtml).

2. Run StrainID

For yeast (sacCer3) samples...

```
cd GenoPipe/StrainID
bash identify-Strain.sh -i /path/to/BAM -o /path/to/output -g /path/to/sacCer3.fa -v /path/to/GenoPipe/StrainID/sacCer3_VCF
```

For human (hg19) samples...

```
cd GenoPipe/StrainID
bash identify-Strain.sh -i /path/to/BAM -o /path/to/output -g /path/to/hg19.fa -v /path/to/GenoPipe/StrainID/hg19_VCF
```

Joe Schmoe Example:

In the following example, GenoPipe, the directory including all the input yeast BAM files, and the new directory for storing DeletionID reports are stored on the Desktop of Joe Schmoe. Filepaths would need to be changed according to a user's preferred directory structure.

```
cd /User/joeschmoe/Desktop/GenoPipe/
cd EpitopeID/utility_scripts/genome_data
bash download_sacCer3_Genome.sh
mv genome.fa /User/joeschmoe/Desktop/GenoPipe/sacCer3.fa
# Run StrainID
cd ../../../StrainID
bash identify-Strain.sh -i /User/joeschmoe/Desktop/mybam -o /User/joeschmoe/Desktop/myreports_SID -g /User/joeschmoe/Desktop/GenoPipe/sacCer3.fa -v /User/joeschmoe/Desktop/GenoPipe/StrainID/sacCer3_VCF
```


Full Joe Schmoe examples' directory structure:

```
/User/joeschmoe/Desktop
  |--GenoPipe
  |   |--EpitopeID
  |   |--DeletionID
  |   |--StrainID
  |--myfastq
  |  |--SampleA_R1.fastq.gz
  |  |--SampleB_R1.fastq.gz
  |  |--SampleB_R2.fastq.gz
	|--mybam
  |  |--SampleA.bam
  |  |--SampleB.bam
  |--myreports_EID
  |  |--SampleA_R1-ID.tab
  |  |--SampleB_R1-ID.tab
  |--myreports_DID
  |  |--SampleA_deletion.tab
  |  |--SampleB_deletion.tab
  |--myreports_SID
     |--SampleA_strain.tab
     |--SampleB_strain.tab
```
