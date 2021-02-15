#!/usr/bin/perl

die "FASTA_File\tOutput_FASTQ.GZ\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;

open(IN, "<$input") or die "Can't open $input for reading!\n";
#open(OUT, ">$output") or die "Can't open $output for writing!\n";
open(OUT, "| gzip -c > $output");

while($line = <IN>) {
    chomp($line);                 # removed indentation
	if((substr $line, 0, 1) eq ">") {
		# line1 write of FASTQ entry
		$ID = substr $line, 1;
		# seems bedtools output adds "()" to the end of the title so you can remove those two chars with the "\n" char
		print OUT "@",(substr $ID, 0, -3),"\n";
		# line2 write
		$line = <IN>;
		chomp($line);
		print OUT $line,"\n";
		# line3 write
		print OUT "+\n";
		# line4 write
		$LENGTH = length($line);
		$PHRED = "";
    	for($x = 0; $x < $LENGTH; $x++) { $PHRED .= "I"; }      # removed indentation, set all quality values to "I", max on Phred+33
		print OUT $PHRED,"\n";
	}
}
close IN;
close OUT;
