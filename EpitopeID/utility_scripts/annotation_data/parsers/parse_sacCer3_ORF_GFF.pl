#! /usr/bin/perl

die "GFF_File\tOutput_GFF\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

# ##gff-version 3
# #date Tue Jan 17 19:50:03 2017
# #
# # Saccharomyces cerevisiae S288C genome (version=R64-2-1)
# #
# # Features from the 16 nuclear chromosomes labeled chrI to chrXVI,
# # plus the mitochondrial genome labeled chrmt.
# #
# # Created by Saccharomyces Genome Database (http://www.yeastgenome.org/)
# #
# # Weekly updates of this file are available for download from:
# # http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff
# #
# # Please send comments and suggestions to sgd-helpdesk@lists.stanford.edu
# #
# # SGD is funded as a National Human Genome Research Institute Biomedical Informatics Resource from
# # the U. S. National Institutes of Health to Stanford University.
# #
# chrI    SGD     chromosome      1       230218  .       .       .       ID=chrI;dbxref=NCBI:;Name=chrI
# chrI    SGD     telomere        1       801     .       -       .       ID=TEL01L;Name=TEL01L;Note=Telomeric%20region%20on%20the%20left%20arm%20of%20Chromosome%20I%3B%20composed%20of%20an%20X%20element%20core%20sequence%2C%20X%20element%20combinatorial%20repeats%2C%20and%20a%20short%20terminal%20stretch%20of%20telomeric%20repeats;display=Telomeric%20region%20on%20the%20left%20arm%20of%20Chromosome%20I;dbxref=SGD:S000028862
# chrI    SGD     X_element       337     801     .       -       .       ID=TEL01L_X_element;Name=TEL01L_X_element;dbxref=SGD:S000028862
# chrI    SGD     X_element_combinatorial_repeat  63      336     .       -       .       ID=TEL01L_X_element_combinatorial_repeat;Name=TEL01L_X_element_combinatorial_repeat;dbxref=SGD:S000028862
# chrI    SGD     telomeric_repeat        1       62      .       -       .       ID=TEL01L_telomeric_repeat;Name=TEL01L_telomeric_repeat;dbxref=SGD:S000028862
# chrI    SGD     gene    335     649     .       +       .       ID=YAL069W;Name=YAL069W;Ontology_term=GO:0003674,GO:0005575,GO:0008150;Note=Dubious%20open%20reading%20frame%3B%20unlikely%20to%20encode%20a%20functional%20protein%2C%20based%20on%20available%20experimental%20and%20comparative%20sequence%20data;display=Dubious%20open%20reading%20frame;dbxref=SGD:S000002143;orf_classification=Dubious
# chrI    SGD     CDS     335     649     .       +       0       Parent=YAL069W_mRNA;Name=YAL069W_CDS;orf_classification=Dubious

$CHR{"chrI"} = "chr1";
$CHR{"chrII"} = "chr2";
$CHR{"chrIII"} = "chr3";
$CHR{"chrIV"} = "chr4";
$CHR{"chrV"} = "chr5";
$CHR{"chrVI"} = "chr6";
$CHR{"chrVII"} = "chr7";
$CHR{"chrVIII"} = "chr8";
$CHR{"chrIX"} = "chr9";
$CHR{"chrX"} = "chr10";
$CHR{"chrXI"} = "chr11";
$CHR{"chrXII"} = "chr12";
$CHR{"chrXIII"} = "chr13";
$CHR{"chrXIV"} = "chr14";
$CHR{"chrXV"} = "chr15";
$CHR{"chrXVI"} = "chr16";
$CHR{"chrmt"} = "chrM";

$line = "";
while($line = <IN>) {
	chomp($line);
	next if($line =~ "#");
	next if($line =~ "Dubious");
	next if($line =~ "Uncharacterized");
	@array = split(/\t/, $line);
	if($array[2] eq "gene") {
		# Overwrite roman numeral to arabic
                $array[0] = $CHR{$array[0]};

		# Create genomic coordinate
		$COORD = $array[0] . ":" . $array[3] . "-" . $array[4];
		
		@ATTR = split(/\;/, $array[8]);
		$NAME = 0;
		for($x = 0; $x <= $#ATTR; $x++) {
			if($ATTR[$x] =~ "gene=") {
				@GENE = split(/\=/, $ATTR[$x]);
	                        @FIRST = split(/\,/, $GENE[1]);
        	                $array[8] = $FIRST[0]  . "|" . $COORD;
				$NAME = 1;
			}
		}
		if($NAME == 0) {
			@NAME = split(/\=/, $ATTR[1]);
			$array[8] = $NAME[1] . "|" . $COORD;
		}

		print OUT join("\t", @array),"\n";

	}
}
close IN;
close OUT;
