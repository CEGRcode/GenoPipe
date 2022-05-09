#! /usr/bin/perl/


die "GENOME_FASTA\tCOORDS\tOutput_GENOME_FASTA\n" unless $#ARGV == 2;
my($genome, $coords, $output) = @ARGV;

# $genome is path to FASTA file of genomic sequence to insert epitope into
# $coords is a string indicating locus to delete (semicolon separated, 1-indexed, remove inclusive interval), e.g. chr2:334386-336818
# $output is path to write new FASTA file to (formatted to use same number of bp per line as $genome FASTA file)

# Works as long as fastaWidth is determined before insertion location is encountered
# Assumes only 1 sequence in input FASTA, loads seq string to be inserted into the genome.

# Determining $coord string value
# GeneCoords Reb1 and Rap1 (1-indexed,inclusive)
#	chr2:334386:336818:-
#	chr14:241689:244172:+
# For Reb1, use 
# 	coord=chr2:334386-336818
# For Rap1, use
#	coord=chr14:241689-244172

# (1) Load in "coord" of sequence to delete from the genome
@array=split /:/, $coords;
$inCHR = $array[0];
@width=split /-/, $array[1];
$inStart = $width[0];
$inEnd = $width[1];

# Check start coordinate has smaller index
if( $inStart > $inEnd ){
	$temp = $inStart;
	$inStart = $inEnd;
	$inEnd = $temp;
}

print "Deletion coordinates loaded: ${inCHR}:${inStart}-${inEnd}\n";


# (2) Copy out genome with sequence from input coordinates removed

open(IN, "<$genome") or die "Can't open $genome for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

$fastaWidth = -999;
$fastaSeq = "";
# $deleteComplete = 0;
$currentChr = "";
$currentBP = 0;

while($line = <IN>) {
	chomp($line);
	# Parse sequence header and reset variables as appropriate
	if((substr $line, 0, 1) eq ">") {
        $currentChr = substr $line, 1;
		if($fastaSeq ne "") { print OUT $fastaSeq,"\n"; $fastaSeq = ""; }
		print OUT $line,"\n";
	# Parse sequence lines of the chromosome where we will insert the input FASTA
	} elsif($currentChr eq $inCHR) {
		# Set fastaWidth for formatting
		if($fastaWidth == -999) { $fastaWidth = length($line); }
		$currentBP += length($line);
		if($inStart > $currentBP) {
			print OUT $line,"\n";   # print until we get to the insertion site block
		} else {
 			if($deleteComplete == 0) {
 				$L = length($line);
 				# Save the sequence before the deletion start site within this block if start occurs within this block
 				$startDiff = $currentBP - $inStart + 1;
 				if($L >= $startDiff){
					$fastaSeq = (substr $line, 0, ( $L - $startDiff ));
 				}
				
				# Save the sequence after the deletion end site within this block if end occurs within this block
				$endDiff = $currentBP - $inEnd;
				if( $inEnd<=$currentBP && $L>$endDiff ){
					$fastaSeq .= (substr $line, ($L - $endDiff), $endDiff);
					$deleteComplete = 1;
				}
				
				# Slice and output a fasta block if fastaSeq long enough
				if( length($fastaSeq)>=$fastaWidth ){
					print OUT (substr $fastaSeq, 0, $fastaWidth),"\n";
	                                $fastaSeq = substr $fastaSeq, $fastaWidth;
				}
			} elsif($deleteComplete == 1) {
				# Append next sequence block and slice out fastaWidth to output
				$fastaSeq .= $line;
				print OUT (substr $fastaSeq, 0, $fastaWidth),"\n";
				$fastaSeq = substr $fastaSeq, $fastaWidth;
 			}
		}
	# Straight Print Sequence lines not having to do with the input chromosome
	} else {
		print OUT $line,"\n";
	}
}
close IN;
close OUT;
