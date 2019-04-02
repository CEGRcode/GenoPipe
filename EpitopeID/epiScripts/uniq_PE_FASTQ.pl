#! /usr/bin/perl

die "FASTA_1\tFASTA_2\n" unless $#ARGV == 1;
my($f1, $f2) = @ARGV;

open(F1, "<$f1") or die "Can't open $f1 for reading!\n";
open(F2, "<$f2") or die "Can't open $f2 for reading!\n";

# >NS500168:386:H37NCBGX9:1:11306:2848:14354 1:N:0:CTGCGCAA+GCGCCTTA
# GCCAGCGACATGGAGGCCCAGAATACCATCCTAGACAGAC
# >NS500168:386:H37NCBGX9:1:11310:5053:16624 1:N:0:CTGCGCAT+GAGCCTTA
# AACGCGTACTATCAAACATTATTGAGTTTTTCGCTTTCAC
# >NS500168:386:H37NCBGX9:1:11310:11477:18526 1:N:0:CTGCGCAT+GAGCCTTA
# AGTCACCGCCGGGTCTCCCGGCCAGCGCCATGGAGGCCCA

%UNIQ = ();
$line1 = "";
$line2 = "";
while($line1 = <F1>) {
	$line2 = <F2>;
	chomp($line1);
	chomp($line2);
	
	$HEADER = $line1 . "|" . $line2;

	$line1 = <F1>;
	$line2 = <F2>;
	chomp($line1);
	chomp($line2);

	$SEQ = $line1 . "|" . $line2;

	$UNIQ{$SEQ} = $HEADER;
}
close F1;
close F2;

foreach $key (keys %UNIQ) {
	@HEAD = split(/\|/, $UNIQ{$key});
	@SEQ = split(/\|/, $key);
	print $HEAD[0],"\n",$SEQ[0],"\n";	
}
