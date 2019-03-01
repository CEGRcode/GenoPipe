#!/bin/bash

# hg38
# This script pulls the latest RefSeq GFF annotation from NCBI and parses it into a format compatible with 
# the epitope ID scripts. Move the resulting files into the appropriate /pwd/gene_FASTA/ when complete

# Required software:
# wget
# bedtools v2.26+
# BWA v0.7.14+

usage()
{
    echo 'get_hg38_transcript_FASTA.sh -f <Genomic FASTA file>'
    echo 'eg: bash get_hg38_transcript_FASTA.sh -f hg38.fa'
    exit
}

if [ "$#" -ne 2 ]
then
    usage
fi

while getopts ":f:" IN; do
    case "${IN}" in
        f)
            FASTA=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${FASTA}" ] ; then
    usage
fi

echo "f = ${FASTA}"

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
gunzip -f GRCh38_latest_genomic.gff.gz
perl parsers/parse_NCBI_RefSeq_Annotation.pl GGRCh38_latest_genomic.gff temp.gff
bedtools getfasta -name -s -fi $FASTA -bed temp.gff -fo gene.fasta
bwa index gene.fasta
rm GGRCh38_latest_genomic.gff temp.gff
