#! /usr/bin/perl

die "GFF_File\tGene_Border_Pad(bp)\tOutput_GFF\n" unless $#ARGV == 2;
my($input, $PAD, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

# W3110   .       .       191     255     .       +       .       thrL
# W3110   .       .       338     2799    .       +       .       thrA
# W3110   .       .       2802    3733    .       +       .       thrB
# W3110   .       .       3735    5020    .       +       .       thrC
# W3110   .       .       5235    5530    .       +       .       yaaX

#chr1    sacCer3.fa      BIN     1501    1556    .       +       .       chr1:1501-1556
#chr1    SGD     geneBorder      1557    1807    .       -       .       PAU8-Cterm|chr1:1557-1807

$line = "";
while($line = <IN>) {
	chomp($line);
	@array = split(/\t/, $line);
	# Overwrite roman numeral to arabic

	# Create genomic coordinate
	$COORD = $array[0] . ":" . $array[3] . "-" . $array[4];
	$GENE = $array[8];
	$array[8] = $GENE . "|" . $COORD;

	# Print out gene info
	print OUT join("\t", @array),"\n";

	# Print out border region of gene given user-specified padding region
	$LOC = "Promoter";
	if($array[6] eq "-") { $LOC = "Cterm"; }
        # Create upstream coordinate
	$UP = $array[3] - $PAD;
	if($UP < 1) { $UP = 1; }
	print OUT "$array[0]\t$array[1]\t$array[2]\t$UP\t$array[3]\t$array[5]\t$array[6]\t$array[7]\t$GENE\-$LOC|$COORD\n";

	$LOC = "Cterm";
        if($array[6] eq "-") { $LOC = "Promoter"; }
	# Create downstream coordinate
	$DOWN = $array[4] + $PAD;
        print OUT "$array[0]\t$array[1]\t$array[2]\t$array[4]\t$DOWN\t$array[5]\t$array[6]\t$array[7]\t$GENE\-$LOC|$COORD\n";
}
close IN;
close OUT;
