#! /usr/bin/perl

die "GFF_File\tOutput_GFF\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

$line = "";
while($line = <IN>) {
	chomp($line);
	@array = split(/\t/, $line);
	$ID = $array[0] . ":" . $array[3] . "-" . $array[4];
	$array[8] = $ID;
	print OUT join("\t", @array),"\n";
}
close IN;
close OUT;
