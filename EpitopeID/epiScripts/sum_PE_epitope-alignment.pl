#! /usr/bin/perl

die "epitope-se.out\talign-pe.out\tOutput_Counts\n" unless $#ARGV == 2;
my($se_epi, $pe_epi, $output) = @ARGV;

open(TAG, "<$se_epi") or die "Can't open $se_epi for reading!\n";

#NS500168:392:H35LFBGX9:1:11106:11140:8946       Extended-Tap
#NS500168:392:H35LFBGX9:1:11202:13231:13791      Extended-Tap
#NS500168:392:H35LFBGX9:1:11205:7905:8901        Extended-Tap
#NS500168:392:H35LFBGX9:1:11207:14438:1246       Extended-Tap
#NS500168:392:H35LFBGX9:1:11209:8225:6951        Extended-Tap
#NS500168:392:H35LFBGX9:1:11210:17172:7664       Extended-Tap

%READ_TAG = ();
$line = "";
while($line = <TAG>) {
	chomp($line);
	@array = split(/\t/, $line);
	$READ_TAG{$array[0]} = $array[1];
}
close TAG;

open(PE, "<$pe_epi") or die "Can't open $pe_epi for reading!\n";

#chr13   255700  255736  NS500168:392:H35LFBGX9:3:12610:3794:5579        60      +       255700  255736  0,0,0   1       36,     0,      chr13   SGD     gene    253848  255800  .       +       .       YAP1|chr13:253848-255800
#chr13   255710  255746  NS500168:392:H35LFBGX9:1:11202:13231:13791      60      +       255710  255746  0,0,0   1       36,     0,      chr13   SGD     gene    253848  255800  .       +       .       YAP1|chr13:253848-255800
#chr13   255649  255685  NS500168:392:H35LFBGX9:2:13211:3874:13665       60      +       255649  255685  0,0,0   1       36,     0,      chr13   SGD     gene    253848  255800  .       +       .       YAP1|chr13:253848-255800
#chr13   255800  255836  NS500168:392:H35LFBGX9:3:21404:17665:2395       60      -       255800  255836  0,0,0   1       36,     0,      chr13   sacCer3.fa      BIN     255801  256001  .       +       .       chr13:255801-256001
#chr13   255838  255874  NS500168:392:H35LFBGX9:1:12308:11416:19001      60      -       255838  255874  0,0,0   1       36,     0,      chr13   sacCer3.fa      BIN     255801  256001  .       +       .       chr13:255801-256001
#chr13   255698  255734  NS500168:392:H35LFBGX9:1:22103:6767:17258       60      +       255698  255734  0,0,0   1       36,     0,      chr13   SGD     gene    253848  255800  .       +       .       YAP1|chr13:253848-255800

%READCOUNT = ();
$line = "";
while($line = <PE>) {
	chomp($line);
	next if((substr $line, 0, 1) eq "@");
	@array = split(/\t/, $line);

	# Set predicted terminus of epitope
	$LOC = "C-term";
	if($array[5] eq "+" && $array[18] eq "-") { $LOC = "N-term"; }
	elsif($array[5] eq "-" && $array[18] eq "+") { $LOC = "N-term"; }
	# Ignore predicted terminus if aligned to a BIN
	if($array[14] eq "BIN" || $array[14] eq "geneBorder") { $LOC = "N/A"; }

	$COUNT = $array[20] . "~" . $READ_TAG{$array[3]} . "~" . $LOC;
	if(exists $READCOUNT{$COUNT} ) { $READCOUNT{$COUNT} = $READCOUNT{$COUNT} + 1; }
	else { $READCOUNT{$COUNT} = 1; }
}
close SAM;

@ARRAY = ();
foreach $key (keys %READCOUNT) { push(@ARRAY, {count => $READCOUNT{$key}, id => $key}); }
@SORT = sort { $$b{'count'} <=> $$a{'count'} } @ARRAY;

open(OUT, ">$output") or die "Can't open $output for writing!\n";
if($#SORT == -1) {
	print OUT "Epitope could not be detected genomically\n";
} else {
	for($x = 0; $x <= $#SORT; $x++) {
		@temparray = split(/\~/, $SORT[$x]{'id'});
		for($y = 0; $y <= $#temparray; $y++) { print OUT "$temparray[$y]\t" }
		print OUT "$SORT[$x]{'count'}\n";;
	}
}
close OUT;
