#! /usr/bin/perl/

die "BED_File\tBORDER_SIZE\tOutput_File\n" unless $#ARGV == 2;
my($input, $BORDER, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for reading!\n";

while($line = <IN>) {
	chomp($line);
	@array = split(/\t/, $line);
	$array[1] -= ($BORDER + 1);
	if($array[1] < 0) { $array[1] = 0; }
	$array[2] += ($BORDER + 1);
	print OUT join("\t", @array),"\n";
}
close IN;
close OUT;

