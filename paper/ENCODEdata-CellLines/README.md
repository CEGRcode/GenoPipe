# StrainID detection of cell lines in ENCODE samples

Run StrainID on ENCODE data and evaluate StrainID's performance.

|  |  |
| :--: | -- |
| Figure 6C  | `ENCODEdata-CellLines/results/ID/` |
| Supplementary Table 10 | `ENCODEdata-CellLines/results/SupplementaryTable10.txt.gz` |

# Reference files

## 210512_sample_metadata.txt
ENCODE metadata was pulled on May 12, 2021 using the `scripts/get_metadata.py` script that pulls all Biosample accessions with classification="cell_line" and whose string matches one of the cell lines we have in our hg19_VCF database. These are used to pull File accessions with type=BAM and assembly=hg19. We did not filter by assay for this analysis. They are saved with all relevant metadata to the `210512_sample_metadata.txt` file.

Command used: `python scripts/get_metadata.py > 210512_sample_metadata.txt`

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
