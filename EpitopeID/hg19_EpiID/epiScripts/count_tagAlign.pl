#! /usr/bin/perl

die "TAG_File\tSAM_File\tOutput_Counts\n" unless $#ARGV == 2;
my($tagfile, $samfile, $output) = @ARGV;

open(TAG, "<$tagfile") or die "Can't open $tagfile for reading!\n";

#NS500168:339:HC53NBGX5:1:11307:23715:7617      Myc-3x
#NS500168:339:HC53NBGX5:1:12202:1249:14252      AID
#NS500168:339:HC53NBGX5:1:12203:26758:14239     AID
#NS500168:339:HC53NBGX5:1:13103:7396:8657       MNase
#NS500168:339:HC53NBGX5:1:13203:15642:2232      AID

%TAG_MATCH = ();
$line = "";
while($line = <TAG>) {
	chomp($line);
	@array = split(/\t/, $line);
	$TAG_MATCH{$array[0]} = $array[1];
	#print $array[1],"\t",$array[0],"\n";
}
close TAG;

open(SAM, "<$samfile") or die "Can't open $samfile for reading!\n";

#@SQ     SN:RAF1 LN:546
#@SQ     SN:REP2 LN:891
#@PG     ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa mem -t 12 sequenceDB/Gene_Sequence/sacCer3_orf_coding.fasta test/temp_tag/tag-reads.fa
#NS500168:339:HC53NBGX5:1:13207:10871:1211       16      SPT3    81      60      40M     *       0       0       AACCACATCACTGATAGAAGATATAGTGAGGGGTCAAGTG        *       NM:i:0  MD:Z:40 AS:i:40 XS:i:0
#NS500168:339:HC53NBGX5:2:13203:22998:11165      16      SPT3    74      60      40M     *       0       0       CCGTAGAAACCACTTCACTGATAGAAGATATAGTGAGGGG        *       NM:i:1  MD:Z:13A26      AS:i:35 XS:i:0
#NS500168:339:HC53NBGX5:2:21108:21301:4158       16      SPT3    74      60      40M     *       0       0       CCGTAGAAACCACATCACTGATAGAAGATATAGTGAGGGG        *       NM:i:0  MD:Z:40 AS:i:40 XS:i:0

%READCOUNT = ();
$line = "";
while($line = <SAM>) {
	chomp($line);
	next if((substr $line, 0, 1) eq "@");
	#print $line,"\n";
	@array = split(/\t/, $line);
	$LOC = "C-term";
	if($array[1] == 16) { $LOC = "N-term"; }
	$COUNT = $array[2] . "~" . $TAG_MATCH{$array[0]} . "~" . $LOC;
	if(exists $READCOUNT{$COUNT} ) { $READCOUNT{$COUNT} = $READCOUNT{$COUNT} + 1; }
	else { $READCOUNT{$COUNT} = 1; }
}
close SAM;

@ARRAY = ();
foreach $key (keys %READCOUNT) { push(@ARRAY, {count => $READCOUNT{$key}, id => $key}); }
@SORT = sort { $$b{'count'} <=> $$a{'count'} } @ARRAY;

open(OUT, ">>$output") or die "Can't open $output for writing!\n";
if($#SORT == -1) {
        print OUT "\nNo epitope detected at genes\n";
} else {
	print OUT "\nGeneID\tEpitopeID\tEpitopeLocation\tEpitopeCount\n";
	for($x = 0; $x <= $#SORT; $x++) {
		@temparray = split(/\~/, $SORT[$x]{'id'});
		for($y = 0; $y <= $#temparray; $y++) { print OUT "$temparray[$y]\t" }
		print OUT "$SORT[$x]{'count'}\n";; 
	}
}
close OUT;
