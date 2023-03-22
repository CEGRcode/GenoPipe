# Downloaded FASTQ files and EpitopeID results go here


## Run EpitopeID on all ENCODE tagged samples


### Update script and download ENCODE samples
```
qsub 00_download_data.pbs
```

### Update script and run EpitopeID on ENCODE samples
```
qsub 01_indexed_runEID.pbs
```

### Compile results into summary report
```
bash 02_tally_results.sh
```


## Make Browser screenshots (Figure 3)

### BrowserData/Genomes
Two synthetic genomes (ID3-eGFP and NR4A1-eGFP) are built by running
```
bash job/03_MakeBrowserData_genomes.sh
#hg19_ID3-Cterm_LAP-tag.fa
#hg19_NR4A1-Cterm_LAP-tag.fa
```

### BrowserData/FASTQ and BrowserData/BAM
Reads filtered to only include read pairs with at least one read mapping to the eGFP tag were aligned to each of the synthetic genomes generated. This was done by running EpitopeID on the samples (with the "clean-up" removal of the directory storing intermeidate files commented out) so that a set of read IDs could be obtained. Then the raw FASTQ files were filtered to only include these read IDs and then aligned.

```
bash job/04_MakeBrowserData_BAM.sh
#hg19_ID3-Cterm_LAP-tag_ENCFF548RTA.bam
#hg19_ID3-Cterm_LAP-tag_ENCFF671VDI.bam
#hg19_NR4A1-Cterm_LAP-tag_ENCFF548RTA.bam
#hg19_NR4A1-Cterm_LAP-tag_ENCFF671VDI.bam

```

### BrowserData/annotations
Annotations of the new ID3-eGFP and NR4A1-eGFP genomes for each of the two ORFs and the LAP-tag are built using the below methodology and saved to `results/BrowserData/annotations/hg19_ID3-Cterm_LAP-tag.bed` and `results/BrowserData/annotations/hg19_NR4A1-Cterm_LAP-tag.bed`.

Gene ORFs were obtained from gencode:

```
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
```

Parse out all start/stop codon info for genes of interest (ID3 and NR4A1).
```
grep '\"ID3\"' ANNOTATIONS/gencode.v19.annotation.gtf_withproteinids |grep 'stop_codon' > ANNOTATIONS/hg19_NoGenotype_features.gtf
grep '\"ID3\"' ANNOTATIONS/gencode.v19.annotation.gtf_withproteinids |grep 'start_codon' >> ANNOTATIONS/hg19_NoGenotype_features.gtf
grep '\"NR4A1\"' ANNOTATIONS/gencode.v19.annotation.gtf_withproteinids |grep 'stop_codon' >> ANNOTATIONS/hg19_NoGenotype_features.gtf
grep '\"NR4A1\"' ANNOTATIONS/gencode.v19.annotation.gtf_withproteinids |grep 'start_codon' >> ANNOTATIONS/hg19_NoGenotype_features.gtf
```

Identify ORF for each gene from gencode entries (select based on IGV hg19 RefSeq annotations).
```
chr1	23885452	23885917	ID3	0	-
chr12	52432493	52452725	NR4A1	0	+
```

Write new BED coordinates shifted as appropriate for the eGFP(LAP) tag
* `hg19_ID3-Cterm_LAP-tag.bed`
* `hg19_NR4A1-Cterm_LAP-tag.bed`


### Screenshots taken

IGV window screenshot coordinate range:

* `hg19_ID3-Cterm_LAP-tag.fa`
  * ID3-locus
    * center -- chr1:(23885452+915)=23886367
    * window -- `chr1:23885367-23887367`
    * modwindow(2kb) -- `chr1:23885001-23886999`
  * NR4A1-locus
    * center -- chr12:52452725
    * window -- `chr12:52451726-52453724`

* `hg19_NR4A1-Cterm_LAP-tag.fa`
  * ID3-locus
    * center -- chr1:23885452
    * window -- `chr1:23884453-23886451`
  * NR4A1-locus
    * center -- chr12:52452725
    * window -- `chr12:52451725-52453725`
    * modwindow(2kb) -- `chr12:52452201-52454199`
