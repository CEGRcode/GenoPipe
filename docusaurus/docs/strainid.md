---
id: strainid
title: StrainID
sidebar_label: StrainID
---

<!-- [strainid-icon]:../static/genopipe-img/strainid-icon.png -->

__StrainID compares a database of VCF files against an aligned BAM file to check for the presence of SNPs in order to determine a likely cell line/strain used in the experiment.__

A report is generated that scores each input BAM with a value for each given VCF file such thataslkdjfal;skjdf alsdkjf

![Figure1C](/genopipe-img/figure1c.png)

## Usage
```bash
identify-Strain.sh -i /path/to/BAM -g /path/to/genome/fASTA -v /path/to/VCF/files -o /path/to/output
```

## Input(`-i`)
StrainID takes a directory pathname for the input and will run StrainID on all the BAM files (`*.bam`) it finds within that directory path when you execute `identify-Strain.sh`.

:::caution
It is extremely important that the BAM files provided match the reference genome build that the VCF files are based off of.
:::

## Reference Genome (`-g`)
A reference genome is required for StrainID to determine the background model of mutation rates. Include the filepath for to the genomic FASTA file for the parameter value. Make sure the genome used matches the genome the VCF coordinates are based off of.

## Reference Files (`-v`)

For StrainID, this is the "database" or directory with all the VCFs used by `identify-Strain.sh`. You will notice that StrainID provides reference files for both yeast (`sacCer3_VCF`) and human (`hg19_VCF`) so you can quickly get started without building up the database from scratch. However, you are free to customize the database by adding your own VCF files.

:::caution
If you add your own VCF files, make sure they match the genome build of both your BAM files and the other VCF files.
:::

Below describes the yeast and human VCF databases and lists which strain files come with the StrainID tool.

* For yeast, the VCF-DB (`sacCer3_VCF/*.vcf`) with strains listed below
* For humans, the VCF-DB (`hg19_VCF/*.vcf`) with strains listed below

|`sacCer3_VCF`|`hg19_VCF`|
|:-----------:|:-------- :|
| BY4741      | A549     |
| BY4742      | HCT116   |
| CEN.PK2-1Ca | HELA    |
| D273-10B    | HepG2     |
| JK9-3d      | K562    |
| RedStar     | LnCap     |
| RM11-1A     | MCF7    |
| SEY6210     | SKnSH     |
| Sigma1278b-10560-6B |  |
| SK1         |   |
| W303        |   |
| Y55         |   |

Whether you use the provided reference files or create your own, the database should use the following directory structure both to ensure that StrainID can find the correct reference files and for simplicity, clarity, and consistency.

```
/name/of/strain/vcf/db
|--strain1.vcf
|--strain2.vcf
|--strain3.vcf
|--strain4.vcf
```


### Customizing your VCFs
Adding strains is very straightforward. Simply add a file in standard [Variant Call Format][vcf-specs] to the `/name/of/strain/vcf/db/` directory.


## Output Report (`-o`)
If StrainID loops through a directory of `n` BAM files, there will be `n` output files written to the directory path indicated by `-o` with file reports named `XXXXXX_strain.tab`.

Below is an example output file where there is a score written to each strain VCF provided with a significant score.
```
  ENCFF000QXV.bam
README.txt	NaN
LnCap.vcf	-4.523561956057013
MCF7.vcf	-3.321928094887362
SKnSH.vcf	2.415037499278844
A549.vcf	-3.0
```

:::warning

fill in output report description

how some scores have lower precision---relatedness between VCFs and

:::

## FAQs

- Q: I used the 14 VCF files in the default sacCer_VCF database but the output only prints out score for 7 of the strains.
  A: First, check the error messages to make sure StrainID isn't quitting part way for any reason (such as hitting resource allocation limits). If StrainID has completed the run with no issues, this can still be the expected behavior. StrainID does not print scores in circumstances such as if there are no reads for alternate allele signal (no)



[yeast-vcf-ref]:https://www.google.com
[human-vcf-ref]:https://www.google.com

[vcf-specs]:https://www.google.com
