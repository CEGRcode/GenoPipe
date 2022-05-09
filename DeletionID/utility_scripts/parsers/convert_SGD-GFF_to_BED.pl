#! /usr/bin/perl/

die "GFF_File\tBED_File\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for reading!\n";

# Set flanking number of nt to trim
#  6 value trims interval to exclude start and stop codons plus an extra codon from each side
$TRIM = 6;

%CHR = ();
$CHR{'chrI'} = "chr1";
$CHR{'chrII'} = "chr2";
$CHR{'chrIII'} = "chr3";
$CHR{'chrIV'} = "chr4";
$CHR{'chrV'} = "chr5";
$CHR{'chrVI'} = "chr6";
$CHR{'chrVII'} = "chr7";
$CHR{'chrVIII'} = "chr8";
$CHR{'chrIX'} = "chr9";
$CHR{'chrX'} = "chr10";
$CHR{'chrXI'} = "chr11";
$CHR{'chrXII'} = "chr12";
$CHR{'chrXIII'} = "chr13";
$CHR{'chrXIV'} = "chr14";
$CHR{'chrXV'} = "chr15";
$CHR{'chrXVI'} = "chr16";

while($line = <IN>) {
	chomp($line);
	next if(substr $line, 0, 1 eq "#");
	@array = split(/\t/, $line);
	if(exists $CHR{$array[0]} && ($array[2] eq "gene" || $array[2] eq "blocked_reading_frame")) {
		@ATTR = split(/\;/, $array[8]);
		@ID = split(/\=/, $ATTR[1]);
		if($array[8] =~ "gene=") { @ID = split(/\=/, $ATTR[2]); }
		# Hardcode skip of ORFs that will never map (closely flank transposon)
		if($ID[1] =~ "YLR157W-D" || $ID[1] =~ "YLR157W-E") { next }
		# Shift start coord for both GFF to BED adjustment and to trim interval size
		$START = $array[3] - 1 + $TRIM;
		# Adjust stop coord to trim interval size
		$STOP = $array[4] - $TRIM;
		
		print OUT "$CHR{$array[0]}\t$START\t$STOP\t$ID[1]\t.\t$array[6]\n";
	}
}
close IN;
close OUT;

