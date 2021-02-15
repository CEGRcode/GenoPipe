#! /usr/bin/perl

die "GENOME_FASTA\tCOORD\tINPUT_FASTA\tOutput_GENOME_FASTA\n" unless $#ARGV == 3;
my($genome, $coord, $input, $output) = @ARGV;

# $genome is path to FASTA file of genomic sequence to insert epitope into
# $coord is a string indicating locus to insert (semicolon separated, 1-indexed, to the right of the given position), e.g. chr2:334386:-
# $input is path to FASTA file with one sequence entry (epitope seq) to insert into genome
# $output is path to write new FASTA file to (formatted to use same number of bp per line as $genome FASTA file)

# Determining $coord string value
# GeneCoords Reb1 and Rap1 (1-indexed,inclusive)
#	chr2:334386:336818:-
#	chr14:241689:244172:+
# For Reb1-Cterm, use 
# 	coord=chr2:334385:-
# 	-1 from gene coord b/c inserted to right of coordinate
# For Rap1-Nterm, use
#	coord=chr14:241688:+
# 	-1 from gene coord b/c inserted to right of coordinate

# Assumes only 1 sequence in input FASTA, loads seq string to be inserted into the genome.


# (1) Load in "input" FASTA sequence to insert into the genome
$SEQ = "";
open(FASTA, "<$input") or die "Can't open $input for reading!\n";
$line = "";
while($line = <FASTA>) {
	chomp($line);
	if(not ((substr $line, 0, 1) eq ">")) {
		$SEQ .= $line;
	}
}
close FASTA;
print "Insertion sequence loaded: ",length($SEQ),"\n";


# (2) Load in "coord" BED coordinate of where to insert into the genome
@coordarray = split /:/,$coord;
$inCHR=$coordarray[0];
$inStart=$coordarray[1];
$inStrand=$coordarray[2];


# (3) Reverse compliment Tag if  inserted onto negative strand
if( $inStrand eq '-' ){
	print "Reverse complement insertion sequence.\n";
	$SEQ =~ s/A/-/g;
	$SEQ =~ s/T/A/g;
	$SEQ =~ s/-/T/g;
	$SEQ =~ s/C/-/g;
	$SEQ =~ s/G/C/g;
	$SEQ =~ s/-/G/g;
	$SEQ = reverse( $SEQ );
}


# (4) Copy out genome with input sequence inserted at designated coordinate
open(IN, "<$genome") or die "Can't open $genome for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";
$fastaWidth = -999;
$fastaSeq = "";
$insertComplete = 0;
$currentChr = "";
$currentBP = 0;
while($line = <IN>) {
	chomp($line);
	# Parse sequence header and reset variables as appropriate
	if((substr $line, 0, 1) eq ">") {
        $currentChr = substr $line, 1;       # Removed indentation
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
			if($insertComplete == 0) {
				# Sandwich the input FASTA between the sequence before and after insertion site within this block
				$fastaSeq = (substr $line, 0, (length($line) - ($currentBP - $inStart)));
				$fastaSeq .= $SEQ;
				$fastaSeq .= (substr $line, (length($line) - ($currentBP - $inStart)));
				# Reformat fastaSeq to fastaWidth sized blocks
				$counter = 0;
				for($x = 0; $x + $fastaWidth < length($fastaSeq); $x += $fastaWidth) {
					print OUT (substr $fastaSeq, $x, $fastaWidth),"\n";
					$counter++;
				}
				# Slice fastaSeq to remainder block after taking out all full line blocks
				$fastaSeq = substr $fastaSeq, ($counter * $fastaWidth);
				$insertComplete = 1;
			} elsif($insertComplete == 1) {
				# Append next sequence block and slice out fastaWidth to output
				$fastaSeq .= $line;
                print OUT (substr $fastaSeq, 0, $fastaWidth),"\n";        # Removed indentation
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
