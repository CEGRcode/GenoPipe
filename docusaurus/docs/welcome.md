---
id: welcome
title: Welcome
sidebar_label: Welcome
---

<img src={require('/../static/genopipe-img/Figure1.png').default} style={{width:60+'%'}}/>

Welcome to the GenoPipe wiki!

GenoPipe is a tool with three modules that check for different aspects of strain background in Next Generation Sequencing(NGS) datasets. Each module tests for a different aspect of strain background: epitope insertions, deletions, and cell line/variant-based background.
  - [EpitopeID][epitopeid-md] identifies and determines the genomic location of epitopes or other inserted sequences relative to genomic loci.
  - [DeletionID][epitopeid-md] identifies significant depletions of aligned NGS tags in the genome relative to a backgroudn model. This tool is set-up to confirm gene knockouts.
  - [StrainID][strainid-md] compares a database of VCF files against an aligned BAM file to check for the presence of SNPs in order to determine a likely cell line/strain used in the experiment.

This is a great QC tool for your data fresh off the sequencer or for checking old published datasets when you can't find a record of the strain background used. These modules work for data from a variety of genomic sequencing assays, not just whole-genome sequencing(WGS) and can be customized to organisms other than yeast and human. Please check the docs and FAQs of the module you're interested in before starting. Confirm that your dataset is appropriate to use with the module.


## Dependencies

* [BWA v0.7.14+][dependency-bwa] (EpitopeID)
* [Samtools v1.7+][dependency-samtools] (EpitopeID)
* [Bedtools v2.26+][dependency-bedtools] (EpitopeID)
* [Python v2.15][dependency-python2] (EpitopeID, DeletionID, StrainID)
  * [Numpy][dependency-numpy] ???
  * [Scipy][dependency-scipy] (DeletionID)
  * [Pysam][dependency-pysam] (StrainID)
* [Perl5][dependency-perl5] (EpitopeID)
* [GNU grep][dependency-gnu-grep] (EpitopeID)
* [wget][dependency-wget] (EpitopeID)

<!-- Epitope List:
* BWA v0.7.14+
* samtools v1.7+
* bedtools v2.26+
* perl5
* python v2.15 with scipy
* GNU grep (BSD grep on MacOSX is >10X slower)

*Epitope get Genome
* wget

Deletion List:
* python v2.15 with scipy

Strain List:
* python v2.15 with pysam
 -->

You can install most of these dependencies as a [conda environment][conda-install] if you have conda setup on your machine using the following command:
```bash
conda create -n my-genopipe-env -c bioconda -c conda-forge bwa scipy samtools bedtools seqkit
```


## Testing


## Set-up: Download & Install

Open your terminal and move to the directory where you want to keep GenoPipe and type the following commands.

```bash
git clone https://github.com/CEGRcode/GenoPipe
cd GenoPipe
```


## Genome Builds
General notes about the defaults used for out-of-the-box usage:

* helper scripts for downloading genome builds described in [EpitopeID guide][epitopeid-md]
  * sacCer3 download + scripts to rename chromosome identifiers to arabic numeral system
  * hg19 download
* Each module's database coord systems are based off the above genome builds (i.e. VCFs in sacCer_VCF and hg19_VCF and BED coordinates in sacCer_Del)



## File Format Standards

* BED Format -- used by DeletionID. See [specifications][bed-specs].
* Variant Call Format (VCF) -- used by StrainID. See [specifications][vcf-specs].



## Getting Help

* Check the specific module guides for a list of FAQs at the bottom. Pay attention to the tips and warnings
* For Bugs: please open an issue on [Github][github-repo] with the following info
    * command you ran with a description of the input files. _Note: we may ask you for a copy of the input files later_
    * entire error message or a description of the problem
    * Software versions and OS you ran the command on



[epitopeid-md]:epitopeid.md
[deletionid-md]:deletionid.md
[strainid-md]:strainid.md


[github-repo]:https://github.com/CEGRcode/GenoPipe
[conda-install]:https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html

[bed-specs]:https://www.google.com
[vcf-specs]:https://www.google.com
[dependency-bwa]:https://www.google.com
[dependency-samtools]:https://www.google.com
[dependency-bedtools]:https://www.google.com
[dependency-perl5]:https://www.google.com
[dependency-python2]:https://www.google.com
[dependency-scipy]:https://www.google.com
[dependency-numpy]:https://www.google.com
[dependency-pysam]:https://www.google.com
[dependency-gnu-grep]:https://www.google.com
[dependency-wget]:https://www.google.com
