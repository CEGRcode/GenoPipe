#! /usr/bin/perl/

die "Genome_FASTA_File\tRead_Length\tResolution (Step)\tOutput_File\n" unless $#ARGV == 3;
my($input, $LENGTH, $STEP, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for reading!\n";

$currentChrom = "";
$currentLine = "";
$currentBP = 0;
while($line = <IN>) {
	chomp($line);
	if($line =~ ">") {
		for($x = 0; $x < length($currentLine) - $LENGTH; $x+=$STEP) {
			$READ = substr $currentLine, $x, $LENGTH;
			print OUT ">$currentChrom\:$x\n",$READ,"\n";
#			print ">$currentChrom\:$x\n",$READ,"\n";
		}
		$currentChrom = substr $line, 1;
                $currentLine = "";
	} else { $currentLine .= $line; }
}
close IN;

for($x = 0; $x < length($currentLine) - $LENGTH; $x+=$STEP) {
	$READ = substr $currentLine, $x, $LENGTH;
	print OUT ">$currentChrom\:$x\n",$READ,"\n";
#        print ">$currentChrom\:$x\n",$READ,"\n";
}
#close OUT;

