#! /usr/bin/perl

die "Length (bp)\tseed(or'-')\tOutput_FASTA_File\n" unless $#ARGV == 2;
my($length, $seed, $output) = @ARGV;

if( $seed ne '-' ){ srand($seed) }

$SEQ = "";
$A = $T = $G = $C = 0;
for($x = 0; $x < $length; $x++) {
# 	sleep(rand);          # eyeball sequences later to make sure there aren't weird polynucleotide stretches
	$nuc = rand;
	if($nuc < 0.25) { $SEQ .= "A"; $A++; }
	elsif($nuc < 0.5) { $SEQ .= "T"; $T++; }
	elsif($nuc < 0.75) { $SEQ .= "G"; $G++; }
	else { $SEQ .= "C"; $C++; }
}

print "A: $A\nT: $T\nG: $G\nC: $C\n";

# Format SEQ variable such that 60bp are displayed per line
$SEQ =~ s/(.{60})/$1\n/gs;

open(OUT, ">$output") or die "Can't open $output for writing!\n";
print OUT ">RANDOM_SEQ\n$SEQ";
close OUT;
