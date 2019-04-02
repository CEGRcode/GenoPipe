#! /usr/bin/perl

die "TAG_File\tOutput_Counts\n" unless $#ARGV == 1;
my($tagfile, $output) = @ARGV;

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
	if(exists $TAG_MATCH{$array[1]}) { $TAG_MATCH{$array[1]} = $TAG_MATCH{$array[1]} + 1; }
	else { $TAG_MATCH{$array[1]} = 1; }
}
close TAG;

@ARRAY = ();
foreach $key (keys %TAG_MATCH) { push(@ARRAY, {count => $TAG_MATCH{$key}, id => $key}); }
@SORT = sort { $$b{'count'} <=> $$a{'count'} } @ARRAY;

open(OUT, ">$output") or die "Can't open $output for writing!\n";
if($#SORT == -1) {
        print OUT "No Tag ID'd\n";
} else {
	print OUT "EpitopeID\tEpitopeCount\n";
	for($x = 0; $x <= $#SORT; $x++) {
		print OUT "$SORT[$x]{'id'}\t$SORT[$x]{'count'}\n";
	}
}
close OUT;
