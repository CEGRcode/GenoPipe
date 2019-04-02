#! /usr/bin/perl

die "sam-header.txt\n" unless $#ARGV == 0;
my($sam) = @ARGV;

open(SAM, "<$sam") or die "Can't open $sam for reading!\n";

# @SQ	SN:chr1	LN:230218
# @SQ	SN:chr2	LN:813184
# @SQ	SN:chr3	LN:316620
# @SQ	SN:chr4	LN:1531933
# @SQ	SN:chr5	LN:576874
# @SQ	SN:chr6	LN:270161
# @SQ	SN:chr7	LN:1090940
# @SQ	SN:chr8	LN:562643
# @SQ	SN:chr9	LN:439888
# @SQ	SN:chr10	LN:745751
# @SQ	SN:chr11	LN:666816
# @SQ	SN:chr12	LN:1078177
# @SQ	SN:chr13	LN:924431
# @SQ	SN:chr14	LN:784333
# @SQ	SN:chr15	LN:1091291
# @SQ	SN:chr16	LN:948066
# @SQ	SN:chrM	LN:85779
# @PG	ID:bwa	PN:bwa	VN:0.7.16a-r1181	CL:bwa mem -t 4 /Users/pughlab/git/GenoPipe/EpitopeID/sacCer3_EpiID/FASTA_genome/genome.fa /Users/pughlab/git/GenoPipe/EpitopeID/20465_Reb1_i5006_BY4741_-_YPD_-_XO_R1/tag-reads.fa

$SIZE = 0;
$line = "";
while($line = <SAM>) {
	chomp($line);
	@array = split(/\s+/, $line);
	next if($array[0] !~ "\@SQ");
	@chromSize = split(/\:/, $array[2]);
	$SIZE += $chromSize[1];
}
close SAM;
print $SIZE,"\n";
