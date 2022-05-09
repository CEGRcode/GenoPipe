#!/bin/bash

# This script makes genomes to simulate from. Two yeast and two human genomes each with a deletion at the specificed loci

DELETE=scripts/delete_interval_from_Genome.pl
[ -d synthetic_genome ] || mkdir synthetic_genome

perl $DELETE ../input/sacCer3.fa chr14:241689-244172 synthetic_genome/sacCer3_Rap1_del.fa
perl $DELETE ../input/sacCer3.fa chr2:334386-336818 synthetic_genome/sacCer3_Reb1_del.fa

perl $DELETE ../input/hg19.fa chr3:184079502-184086383 synthetic_genome/hg19_CTCF_del.fa
perl $DELETE ../input/hg19.fa chr16:67596310-67673088 synthetic_genome/hg19_POLR2H_del.fa
