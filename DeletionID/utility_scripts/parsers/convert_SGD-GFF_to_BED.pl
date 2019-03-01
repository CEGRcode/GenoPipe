#! /usr/bin/perl/

die "GFF_File\tBED_File\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for reading!\n";

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
	if(exists $CHR{$array[0]} && $array[2] eq "gene" && $array[8] =~ "orf_classification=Verified") {
		#Parent=YBR050C_mRNA;Name=YBR050C_CDS;orf_classification=Verified
		@ATTR = split(/\;/, $array[8]);
		@ID = split(/\=/, $ATTR[1]);
		if($array[8] =~ "gene=") { @ID = split(/\=/, $ATTR[2]); }
		$START = $array[3] - 1;
		print OUT "$CHR{$array[0]}\t$START\t$array[4]\t$ID[1]\t.\t$array[6]\n";
	}
}
close IN;
close OUT;

