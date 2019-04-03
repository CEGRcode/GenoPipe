#! /usr/bin/perl

die "GFF.GZ_File\tOutput_BED\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;
open(IN, "gunzip -c $input |") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

# ##gff-version 3
# #!gff-spec-version 1.21
# #!processor NCBI annotwriter
# #!genome-build GRCh37.p13
# #!genome-build-accession NCBI_Assembly:GCF_000001405.25
# #!annotation-date 
# #!annotation-source 
# ##sequence-region NC_000001.10 1 249250621
# ##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606
# NC_000001.10	RefSeq	region	1	249250621	.	+	.	ID=id0;Dbxref=taxon:9606;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA
# NC_000001.10	BestRefSeq	gene	11874	14409	.	+	.	ID=gene0;Dbxref=GeneID:100287102,HGNC:HGNC:37102;Name=DDX11L1;description=DEAD/H-box helicase 11 like 1;gbkey=Gene;gene=DDX11L1;gene_biotype=misc_RNA;pseudo=true
# NC_000001.10	BestRefSeq	transcript	11874	14409	.	+	.	ID=rna0;Parent=gene0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;Name=NR_046018.2;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
# NC_000001.10	BestRefSeq	exon	11874	12227	.	+	.	ID=id1;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
# NC_000001.10	BestRefSeq	exon	12613	12721	.	+	.	ID=id2;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
# NC_000001.10	BestRefSeq	exon	13221	14409	.	+	.	ID=id3;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2
# NC_000001.10	BestRefSeq	gene	14362	29370	.	-	.	ID=gene1;Dbxref=GeneID:653635,HGNC:HGNC:38034;Name=WASH7P;description=WAS protein family homolog 7 pseudogene;gbkey=Gene;gene=WASH7P;gene_biotype=misc_RNA;gene_synonym=FAM39F,WASH5P;pseudo=true
# NC_000001.10	BestRefSeq	transcript	14362	29370	.	-	.	ID=rna1;Parent=gene1;Dbxref=GeneID:653635,Genbank:NR_024540.1,HGNC:HGNC:38034;Name=NR_024540.1;gbkey=misc_RNA;gene=WASH7P;product=WAS protein family homolog 7 pseudogene;transcript_id=NR_024540.1
# NC_000001.10	BestRefSeq	exon	29321	29370	.	-	.	ID=id4;Parent=rna1;Dbxref=GeneID:653635,Genbank:NR_024540.1,HGNC:HGNC:38034;gbkey=misc_RNA;gene=WASH7P;product=WAS protein family homolog 7 pseudogene;transcript_id=NR_024540.1
# NC_000001.10	BestRefSeq	exon	24738	24891	.	-	.	ID=id5;Parent=rna1;Dbxref=GeneID:653635,Genbank:NR_024540.1,HGNC:HGNC:38034;gbkey=misc_RNA;gene=WASH7P;product=WAS protein family homolog 7 pseudogene;transcript_id=NR_024540.1

$CHR{"NC_000067.6"} = "chr1";
$CHR{"NC_000068.7"} = "chr2";
$CHR{"NC_000069.6"} = "chr3";
$CHR{"NC_000070.6"} = "chr4";
$CHR{"NC_000071.6"} = "chr5";
$CHR{"NC_000072.6"} = "chr6";
$CHR{"NC_000073.6"} = "chr7";
$CHR{"NC_000074.6"} = "chr8";
$CHR{"NC_000075.6"} = "chr9";
$CHR{"NC_000076.6"} = "chr10";
$CHR{"NC_000077.6"} = "chr11";
$CHR{"NC_000078.6"} = "chr12";
$CHR{"NC_000079.6"} = "chr13";
$CHR{"NC_000080.6"} = "chr14";
$CHR{"NC_000081.6"} = "chr15";
$CHR{"NC_000082.6"} = "chr16";
$CHR{"NC_000083.6"} = "chr17";
$CHR{"NC_000084.6"} = "chr18";
$CHR{"NC_000085.6"} = "chr19";
$CHR{"NC_000086.7"} = "chrX";
$CHR{"NC_000087.7"} = "chrY";
$CHR{"NC_005089.1"} = "chrM";

$line = "";
while($line = <IN>) {
	chomp($line);
	next if($line =~ "#");
	@array = split(/\t/, $line);
	if($line =~ "gene_biotype=protein_coding" && exists $CHR{$array[0]}) {
		# Overwrite NCBI formal chromosome name with common name
                $array[0] = $CHR{$array[0]};
	
		$START = $array[3];
		$STOP = $array[4] - 1;

		# Create genomic coordinate
		$COORD = $array[0] . ":" . $array[3] . "-" . $array[4];
		
		@ATTR = split(/\;/, $array[8]);
		$NAME = 0;
		for($x = 0; $x <= $#ATTR; $x++) {
			if($ATTR[$x] =~ "Name=") {
				@GENE = split(/\=/, $ATTR[$x]);
	                        @FIRST = split(/\,/, $GENE[1]);
        	                $array[8] = $FIRST[0]  . "|" . $COORD;
				$NAME = 1;
			}
		}
		
		print OUT "$array[0]\t$START\t$STOP\t$array[8]\t.\t$array[6]\n";
#		print OUT join("\t", @array),"\n";
	}
}
close IN;
close OUT;
