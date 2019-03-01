# GenoPipe

## Toolkit for characterizing the genotype of NGS datasets

There are 3 primary modules for genotype identification:

### EpitopeID

Identify and determine the genomic location of epitopes relative to genomic loci

### DeletionID

Identify signficant depletion of aligned NGS tags in the genome relative to a background model. This module is useful for confirming gene knockouts.

### StrainID

Compare a database of VCF files against an aligned BAM file to check for the presence of SNPs in order to determine a likely cell line/strain used in the experiment
