#! /usr/bin/perl

die "GFF_File\tOutput_File\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;

# #!gff-spec-version 1.21
# #!processor NCBI annotwriter
# #!genome-build GRCh37.p13
# #!genome-build-accession NCBI_Assembly:GCF_000001405.25
# #!annotation-date 
# #!annotation-source 
# ##sequence-region NC_000001.10 1 249250621
# ##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606
# NC_000001.10    RefSeq  region  1       249250621       .       +       .       ID=id0;Dbxref=taxon:9606;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA
# NC_000001.10    BestRefSeq      gene    11874   14409   .       +       .       ID=gene0;Dbxref=GeneID:100287102,HGNC:HGNC:37102;Name=DDX11L1;description=DEAD/H-box helicase 11 like 1;gbkey=Gene;gene=DDX11L1;gene_biotype=misc_RNA;pseudo=true
# NC_000001.10    BestRefSeq      transcript      11874   14409   .       +       .       ID=rna0;Parent=gene0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;Name=NR_046018.2;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
# NC_000001.10    BestRefSeq      exon    11874   12227   .       +       .       ID=id1;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
# NC_000001.10    BestRefSeq      exon    12613   12721   .       +       .       ID=id2;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
# NC_000001.10    BestRefSeq      exon    13221   14409   .       +       .       ID=id3;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2

%CHR = ();
$CHR{"NC_000001.10"} = "chr1";
$CHR{"NC_000002.11"} = "chr2";
$CHR{"NC_000003.11"} = "chr3";
$CHR{"NC_000004.11"} = "chr4";
$CHR{"NC_000005.9"} = "chr5";
$CHR{"NC_000006.11"} = "chr6";
$CHR{"NC_000007.13"} = "chr7";
$CHR{"NC_000008.10"} = "chr8";
$CHR{"NC_000009.11"} = "chr9";
$CHR{"NC_000010.10"} = "chr10";
$CHR{"NC_000011.9"} = "chr11";
$CHR{"NC_000012.11"} = "chr12";
$CHR{"NC_000013.10"} = "chr13";
$CHR{"NC_000014.8"} = "chr14";
$CHR{"NC_000015.9"} = "chr15";
$CHR{"NC_000016.9"} = "chr16";
$CHR{"NC_000017.10"} = "chr17";
$CHR{"NC_000018.9"} = "chr18";
$CHR{"NC_000019.9"} = "chr19";
$CHR{"NC_000020.10"} = "chr20";
$CHR{"NC_000021.8"} = "chr21";
$CHR{"NC_000022.10"} = "chr22";
$CHR{"NC_000023.10"} = "chrX";
$CHR{"NC_000024.9"} = "chrY";

open(GFF, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

$line = "";
while($line = <GFF>) {
	chomp($line);
	next if(substr $line, 0, 1 eq "#");

	@array = split(/\t/, $line);
	if($array[2] eq "gene") {
		@ATTR = split(/\;/, $array[8]);
		@TYPE = split(/\=/, $ATTR[6]);
		if($TYPE[1] eq "protein_coding" && exists $CHR{$array[0]}) {
			@GENEID = split(/\=/, $ATTR[5]);
			print OUT $CHR{$array[0]},"\t",$array[1],"\t",$GENEID[1];
			for($x = 3; $x <= $#array; $x++) { print OUT "\t",$array[$x]; }
			print OUT "\n";
#			print OUT $line,"\t",$GENEID[1],"\n";
		}
	}
}
close GFF;
close OUT;
