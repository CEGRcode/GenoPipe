#! /usr/bin/perl

die "Genomic_FASTA\t#_BED\tseed\tOutput_Name\n" unless $#ARGV == 3;
my($genomic,$NUM,$seed,$output) = @ARGV;
open(OUT1, "| gzip -c > $output\_R1.bed.gz");
open(OUT2, "| gzip -c > $output\_R2.bed.gz");
#open(OUT1, ">$output\_R1.bed") or die "Can't open $output\_R1.bed for writing!\n";
#open(OUT2, ">$output\_R2.bed") or die "Can't open $output\_R1.bed for writing!\n";

if( $seed ne '-' ){ srand($seed) }

$genomeSize = 0;
%CHR;
$READLENGTH = 40;
$PEINSERT = 200; # int(rand($PEINSERT)) + 100 creates random insert from 100-300 bp

# Get chr sequence lengths and genome sizes
open(GEN, "<$genomic") or die "Can't open $genomic for reading!\n";
$line = "";
$currentLength = 0;
$currentChr = "";
while($line = <GEN>) {
	chomp($line);
	if((substr $line, 0, 1) eq ">") {
		if($currentChr ne "") {	$CHR{$currentChr} = $currentLength; }
		$currentChr = substr $line, 1;
		$currentLength = 0;
	} else {
		$currentLength += length($line);
		$genomeSize += length($line);
	}
}
$CHR{$currentChr} = $currentLength;
close GEN;

# Simulate reads with length set above(40bp) and 100 to 350 bp insert size
# R1 and R2 coord info stored in read ID of each read in the output R1 and R2 files.
print "Genomesize: ",$genomeSize,"\n";
for($x = 0; $x < $NUM; $x++) {
	$COORD = int(rand($genomeSize));
	for $key (keys %CHR) {
		if($COORD - $CHR{$key} > 0) { $COORD -= $CHR{$key}; }
		else {
			$READID = $CHR{$key} . ":" . $COORD;       # ASK WILL: is this line necessary? variable overwritten a couple lines down and unclear if this var is used before then
			if(rand(1) < 0.5) { $DIR = "+"; }
			else { $DIR = "-"; }
			$R1_END = ($COORD + $READLENGTH);
			$READID = $key . ":" . $COORD . "-" . $R1_END . "," . $DIR;
			
			$COORD2 = $COORD + int(rand(250) + 100);     # ASK WILL: did he want to substitute 250bp with $PEINSERT variable specified above?
			$R2_END = ($COORD2 + $READLENGTH);
			$DIR2 = "+";
			if($DIR eq "+") { $DIR2 = "-"; }
            
            $READID2 = $key . ":" . $COORD2 . "-" . $R2_END . "," . $DIR2;        # Removed indentation here
                        
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
