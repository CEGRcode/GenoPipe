#! /usr/bin/perl
#
#==> Reb1_10K_1000_R1-ID.tab <==
#EpitopeID       EpitopeCount
#RANDOM_SEQ      1
#
#GeneID  EpitopeID       EpitopeLocation EpitopeCount    pVal
#REB1-Cterm|chr2:334136-334386           N/A     1       4.112736692552819e-05
#
#==> Reb1_10K_100_R1-ID.tab <==
#EpitopeID       EpitopeCount
#RANDOM_SEQ      2
#
#GeneID  EpitopeID       EpitopeLocation EpitopeCount    pVal
#Epitope could not be detected genomically
#
#==> Reb1_10K_101_R1-ID.tab <==
#EpitopeID       EpitopeCount
#RANDOM_SEQ      2


die "Stats_File\tGeneName\tOutput_Stats\n" unless $#ARGV == 2;
my($input, $correct, $output) = @ARGV;

# $input	use head with wildcards to get contents of all ID results
# $correct	string that indicates correct call (gene name or tile name of insert locus)
# $output	file of stats summary on simulations

$CORRECT = uc $correct;

open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

$line = "";
%EPITOPE = ();
%CORRECTCALL = ();

while($line = <IN>) {
	chomp($line);
	@array = split(/\s+/, $line);
	if($array[0] eq "EpitopeID") {
		$line = <IN>;
		chomp($line);
		@array = split(/\s+/, $line);
		if($array[0] eq "RANDOM_SEQ") {
			$EPITOPE{$array[1]} = $EPITOPE{$array[1]} + 1;
		} else {
			$EPITOPE{"0"} = $EPITOPE{"0"} + 1;
		}
	}
	if($array[0] eq "GeneID") {
		$CALL = 0;
		$line = <IN>;
		chomp($line);
		while($line =~ $CORRECT) {
			$CALL++;
			$line = <IN>;
	                chomp($line);
		}
		$CORRECTCALL{$CALL} = $CORRECTCALL{$CALL} + 1;
	}
}

close IN;

print OUT "EpitopeID\tFrequency\n";
for $key (keys %EPITOPE) {
	print OUT $key,"\t",$EPITOPE{$key},"\n";

}

print OUT "\nCorrectCall\tFrequency\n";
for $key (keys %CORRECTCALL) {
        print OUT $key,"\t",$CORRECTCALL{$key},"\n";

}

close OUT;
