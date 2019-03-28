#! /usr/bin/perl

die "ALIGN-PE_File\tOutput_File\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;

#chr13   255700  255736  NS500168:392:H35LFBGX9:3:12610:3794:5579        60      +       255700  255736  0,0,0   1       36,     0,      chr13   SGD     gene    253848  255800  .       +       .       YAP1|chr13:253848-255800
#chr13   255710  255746  NS500168:392:H35LFBGX9:1:11202:13231:13791      60      +       255710  255746  0,0,0   1       36,     0,      chr13   SGD     gene    253848  255800  .       +       .       YAP1|chr13:253848-255800
#chr13   255649  255685  NS500168:392:H35LFBGX9:2:13211:3874:13665       60      +       255649  255685  0,0,0   1       36,     0,      chr13   SGD     gene    253848  255800  .       +       .       YAP1|chr13:253848-255800
#chr13   255800  255836  NS500168:392:H35LFBGX9:3:21404:17665:2395       60      -       255800  255836  0,0,0   1       36,     0,      chr13   sacCer3.fa      BIN     255801  256001  .       +       .       chr13:255801-256001
#chr13   255838  255874  NS500168:392:H35LFBGX9:1:12308:11416:19001      60      -       255838  255874  0,0,0   1       36,     0,      chr13   sacCer3.fa      BIN     255801  256001  .       +       .       chr13:255801-256001
#chr13   255698  255734  NS500168:392:H35LFBGX9:1:22103:6767:17258       60      +       255698  255734  0,0,0   1       36,     0,      chr13   SGD     gene    253848  255800  .       +       .       YAP1|chr13:253848-255800
#chr13   255865  255901  NS500168:392:H35LFBGX9:2:21110:3765:13661       60      -       255865  255901  0,0,0   1       36,     0,      chr13   sacCer3.fa      BIN     255801  256001  .       +       .       chr13:255801-256001

open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for reading!\n";

$line = "";
while($line = <IN>) {
	chomp($line);
	@array = split(/\t/, $line);

	$FIVEPRIME = $array[6];
	if($array[5] eq "-") { $FIVEPRIME = $array[7]; }
	if($FIVEPRIME >= $array[15] && $FIVEPRIME <= $array[16]) { print OUT $line,"\n"; }
}
close SAM;

