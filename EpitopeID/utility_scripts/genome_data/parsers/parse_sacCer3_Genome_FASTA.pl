#! /usr/bin/perl

die "FASTA_File\tOutput_FASTA\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

# >ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]
# CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACA
# CATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTT
# ACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC

$CHROM{">ref|NC_001133|"} = "chr1";
$CHROM{">ref|NC_001134|"} = "chr2";
$CHROM{">ref|NC_001135|"} = "chr3";
$CHROM{">ref|NC_001136|"} = "chr4";
$CHROM{">ref|NC_001137|"} = "chr5";
$CHROM{">ref|NC_001138|"} = "chr6";
$CHROM{">ref|NC_001139|"} = "chr7";
$CHROM{">ref|NC_001140|"} = "chr8";
$CHROM{">ref|NC_001141|"} = "chr9";
$CHROM{">ref|NC_001142|"} = "chr10";
$CHROM{">ref|NC_001143|"} = "chr11";
$CHROM{">ref|NC_001144|"} = "chr12";
$CHROM{">ref|NC_001145|"} = "chr13";
$CHROM{">ref|NC_001146|"} = "chr14";
$CHROM{">ref|NC_001147|"} = "chr15";
$CHROM{">ref|NC_001148|"} = "chr16";
$CHROM{">ref|NC_001224|"} = "chrM";

$line = "";
while($line = <IN>) {
	chomp($line);
	@array = split(/\s+/, $line);
	if($array[0] =~ ">") {
		print OUT ">$CHROM{$array[0]}\n";
	} else { print OUT $line,"\n"; }
}
close IN;
close OUT;
