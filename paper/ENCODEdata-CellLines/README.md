# StrainID detection of cell lines in ENCODE samples

Run StrainID on ENCODE data and evaluate StrainID's performance.

|  |  |
| :--: | -- |
| Figure 6C  | `ENCODEdata-CellLines/results/ID/` |
| Supplementary Figure 2  | `ENCODEdata-CellLines/results/hg38_ID/` |
| Supplementary Table 10 | `ENCODEdata-CellLines/results/SupplementaryTable10.txt.gz` |
| Supplementary Table 11 | `ENCODEdata-CellLines/results/SupplementaryTable11.txt.gz` |

# Reference files

## 210512_sample_metadata.txt
ENCODE metadata was pulled on May 12, 2021 using the `scripts/get_metadata.py` script that pulls all Biosample accessions with classification="cell_line" and whose string matches one of the cell lines we have in our hg19_VCF database. These are used to pull File accessions with type=BAM and assembly=hg19. We did not filter by assay for this analysis. They are saved with all relevant metadata to the `210512_sample_metadata.txt` file.

Command used: `python scripts/get_metadata.py > 210512_sample_metadata.txt`

## 230612_hg38_100-TF-ChIP-seq.txt
ENCODE metadata was pulled on June 12, 2023 according to the retrieval URL in the header. The filter criteria ensured files were BAM formats using the hg38 genome build, from the TF ChIP-seq assay, and from a cell line background in the core `hg19_VCF` set. Only the top 100 BAM files were saved here.

# Download ENCODE CellLine data and run through StrainID
Use the 210512_sample_metadata.txt file with ENCODE accessions to download and process the
data to the `results` directory. Then run the data through StrainID.

## Download BAM files
```
qsub job/00_download_data.pbs
```

## Filter BAM files

```
qsub job/01_filter_and_count.pbs
```

## Run StrainID
```
qsub job/02_indexed_runSID.pbs
```

## Compile the results with the metadata information
Evaluate the accuracy of StrainID on real data by merging the metadata with the StrainID results.
```
python scripts/analyze_encode_results.py -i results/ID -m 210512_sample_metadata.txt -o results/SupplementaryTable10.txt
python scripts/build_violinscatter.py -i results/SupplementaryTable10.txt -o results/Figure6C.png
```

You can even filter the visualization to show the results for specific assays by using the `-a` option flag
```
python scripts/build_violinscatter.py -i results/SupplementaryTable10.txt -o results/Figure6C_ChIP-seq.png -a "ChIP-seq"
python scripts/build_violinscatter.py -i results/SupplementaryTable10.txt -o results/Figure6C_CAGE.png -a "CAGE"
python scripts/build_violinscatter.py -i results/SupplementaryTable10.txt -o results/Figure6C_small-RNA-seq.png -a "small RNA-seq"
```

# Characterize StrainID scores from 100 ChIP-seq datasets
To better understand the distribution of scores output by StrainID for "correct" and "incorrect" strain backgrounds, we ran StrainID on our largest reference set of over 1,000 cell lines (`hg38_DepMap`) in a sample of TF-ChIP-seq datasets from ENZCODE.

## Download data & Run StrainID
```
qsub job/02_indexed_runSID.pbs
```

## Compile the results with the metadata information
Compile the results into a table and plot the scores as two overlapping histograms.
```
python scripts/merge_sidscores.py -i results/hg38_ID/ -m 230612_hg38_100-TF-ChIP-seq.txt -o results/SupplementaryTable11.txt
python scripts/build_sidscore_histogram.py -i results/hg38_ID/ -o results/SupplementaryFigure2.png
```
