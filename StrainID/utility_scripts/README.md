
### Model.csv
Downloaded from https://depmap.org/portal/download/all/?releasename=DepMap+Public+23Q2&filename=Model.csv on May 14, 2023


### generate_hg38_VariantDB.sh
Download all DepMap (Broad) variants based on `Model.csv` accessions and cell line names.

### generate_sacCer3_VariantDB.sh
Download all Song et al (2015) variants and perform BedTools intersect for unique substitutions (better uniquely characterize strain) and format into VCF files for StrainID.
