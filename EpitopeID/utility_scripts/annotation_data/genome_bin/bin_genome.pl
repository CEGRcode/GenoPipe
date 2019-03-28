#! /usr/bin/perl

die "Genome_FASTA_File\tBin_Resolution\tOutput_GFF\n" unless $#ARGV == 2;
my($input, $RES, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

@INPUT = split(/\//, $input);

$line = "";
$CURRENTCHROM = "";
$LEN = 0;

while($line = <IN>) {
	chomp($line);
	if($line =~ ">") {
		if($LEN > 0) {
			for($START = 1; $START <= $LEN; $START += $RES) {
				$STOP = $START + $RES;
				if($STOP > $LEN) { $STOP = $LEN; }
				print OUT "$CURRENTCHROM\t$INPUT[$#INPUT]\tBIN\t$START\t$STOP\t.\t+\t.\t$CURRENTCHROM\:$START\-$STOP\n";
			}
		}
		$LEN = 0;
		$CURRENTCHROM = substr $line, 1;
	} else { $LEN += length($line);	}
}
close IN;

if($LEN > 0) {
	for($START = 1; $START <= $LEN; $START += $RES) {
	        $STOP = $START + $RES;
		if($STOP > $LEN) { $STOP = $LEN; }
                print OUT "$CURRENTCHROM\t$INPUT[$#INPUT]\tBIN\t$START\t$STOP\t.\t+\t.\t$CURRENTCHROM\:$START\-$STOP\n";
	}
}
close OUT;
