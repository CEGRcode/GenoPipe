#! /usr/bin/perl

die "Genome_FAI\t#_BED\tseed\tOutput_Name\n" unless $#ARGV == 3;
my($genomic,$NUM,$seed,$output) = @ARGV;
open(OUT1, "| gzip -c > $output\_R1.bed.gz");
open(OUT2, "| gzip -c > $output\_R2.bed.gz");
#open(OUT1, ">$output\_R1.bed") or die "Can't open $output\_R1.bed for writing!\n";
#open(OUT2, ">$output\_R2.bed") or die "Can't open $output\_R1.bed for writing!\n";

if( $seed ne '-' ){ srand($seed) }

$genomeSize = 0;
%CHR;
$READLENGTH = 40;
$PEINSERT = 250;

# Get chr sequence lengths and genome sizes
open(GEN, "<$genomic") or die "Can't open $genomic for reading!\n";
$line = "";
$currentLength = 0;
$currentChr = "";
while($line = <GEN>) {
	chomp($line);
	my @array = split("\t", $line);
	$CHR{$array[0]} = $array[1];
	$genomeSize += $array[1];
}
close GEN;


# Simulate reads with length set above(40bp) and 100 to 350 bp insert size
# R1 and R2 coord info stored in read ID of each read in the output R1 and R2 files.
print "Genomesize: ",$genomeSize,"\n";
for($x = 0; $x < $NUM; $x++) {
  # Set fragment/R1 starting coordinate
	$COORD = int(rand($genomeSize));  # sample random coordinate
	for $key (keys %CHR) {
		if($COORD - $CHR{$key} > 0) { $COORD -= $CHR{$key}; }
		else {
			# Set strand
			if(rand(1) < 0.5) { $DIR = "+"; }   # sample random strand
			else { $DIR = "-"; }
			# Define Read1 interval
			$R1_END = ($COORD + $READLENGTH);
			$READID = $key . ":" . $COORD . "-" . $R1_END . "," . $DIR;

			# Set R2 strand
			$DIR2 = "+";
			if($DIR eq "+") { $DIR2 = "-"; }
			# Define Read2 interval
			$COORD2 = $COORD + int(rand($PEINSERT) + 100);    # sample random insert size
			$R2_END = ($COORD2 + $READLENGTH);
			$READID2 = $key . ":" . $COORD2 . "-" . $R2_END . "," . $DIR2;

			# Unique read id shared between R1 and R2
			$COMBINED = $READID . "|" . $READID2;					# ASK WILL: possibly add x variable to the unique identifier here? (unlikely both R1 and R2 both match...could be considered duplicates)

			# Write read info to BED format (R1 written to OUT1 and R2 written to OUT2)
			if($COORD > 0 && $R2_END < $CHR{$key}) {
				print OUT1 "$key\t$COORD\t$R1_END\t$COMBINED\t.\t$DIR\n";         # Removed indentation here
				print OUT2 "$key\t$COORD2\t$R2_END\t$COMBINED\t.\t$DIR2\n";
			} else { $x--; }   # skip and redo in another simulation if this runs off the end of the chromosome
			last;
		}
	}
#	print $COORD,"\n";
}

close OUT1;
close OUT2;
