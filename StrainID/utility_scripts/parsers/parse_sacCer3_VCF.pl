#! /usr/bin/perl

die "VCF_File\tOutput_VCF\n" unless $#ARGV == 1;
my($input, $output) = @ARGV;
open(IN, "<$input") or die "Can't open $input for reading!\n";
open(OUT, ">$output") or die "Can't open $output for writing!\n";

# ##fileformat=VCFv4.1
# ##FILTER=<ID=LowQual,Description="Low quality">
# ##FILTER=<ID=SnpCluster,Description="SNPs found in clusters">
# ##FILTER=<ID=filter,Description="(AB ?: 0) > 0.75 || QUAL < 50.0 || DP > 360 || SB > -0.1 || MQ0>=4">
# ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
# ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
# ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
# ##GATKCommandLine=<ID=VariantFiltration,CommandLineOptions="analysis_type=VariantFiltration input_file=[] read_buffer_size=null phone_home=AWS gatk_key=null tag=NA read_filter=[] intervals=null excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=/srv/gs1/projects/cherry/sacCer3/sacCer3.fa nonDeterministicRandomSeed=false disableDithering=false maxRuntime=-1 maxRuntimeUnits=MINUTES downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=1000 baq=OFF baqGapOpenPenalty=40.0 fix_misencoded_quality_scores=false allow_potentially_misencoded_quality_scores=false useOriginalQualities=false defaultBaseQualities=-1 performanceLog=null BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 globalQScorePrior=-1.0 allow_bqsr_on_reduced_bams_despite_repeated_warnings=false validation_strictness=SILENT remove_program_records=false keep_program_records=false sample_rename_mapping_file=null unsafe=null disable_auto_index_creation_and_locking_when_reading_rods=false num_threads=1 num_cpu_threads_per_data_thread=1 num_io_threads=0 monitorThreadEfficiency=false num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false variant_index_type=DYNAMIC_SEEK variant_index_parameter=-1 logging_level=INFO log_to_file=null help=false version=false variant=(RodBinding name=variant source=/srv/gs1/projects/cherry/giltae/strains/hugeseq/BY4741/chrI.snp.gatk.vcf) mask=(RodBinding name= source=UNBOUND) out=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub no_cmdline_in_header=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub bcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub filterExpression=[(AB ?: 0) > 0.75 || QUAL < 50.0 || DP > 360 || SB > -0.1 || MQ0>=4] filterName=[filter] genotypeFilterExpression=[] genotypeFilterName=[] clusterSize=3 clusterWindowSize=10 maskExtension=0 maskName=Mask filterNotInMask=false missingValuesInExpressionsShouldEvaluateAsFailing=false invalidatePreviousFilters=false filter_reads_with_N_cigar=false filter_mismatching_base_and_quals=false filter_bases_not_stored=false",Date="Tue Jan 14 23:14:32 PST 2014",Epoch=1389770072639,Version=2.8-1-g932cd3a>
# ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
# ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
# ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
# ##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
# ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
# ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
# ##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
# ##INFO=<ID=Dels,Number=1,Type=Float,Description="Fraction of Reads Containing Spanning Deletions">
# ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
# ##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
# ##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
# ##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
# ##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
# ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
# ##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
# ##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
# ##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
# ##INFO=<ID=RPA,Number=.,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
# ##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
# ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
# ##INFO=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
# ##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
# ##UnifiedGenotyper="analysis_type=UnifiedGenotyper input_file=[/srv/gs1/projects/cherry/giltae/strains/hugeseq/BY4741/chrI.recal.bam] read_buffer_size=null phone_home=STANDARD gatk_key=null read_filter=[] intervals=null excludeIntervals=null interval_set_rule=UNION interval_merging=ALL interval_padding=0 reference_sequence=/srv/gs1/projects/cherry/sacCer3/sacCer3.fa nonDeterministicRandomSeed=false downsampling_type=BY_SAMPLE downsample_to_fraction=null downsample_to_coverage=1000 baq=CALCULATE_AS_NECESSARY baqGapOpenPenalty=40.0 performanceLog=null useOriginalQualities=false BQSR=null quantize_quals=0 disable_indel_quals=false emit_original_quals=false preserve_qscores_less_than=6 defaultBaseQualities=-1 validation_strictness=LENIENT remove_program_records=false keep_program_records=false unsafe=null num_threads=1 num_cpu_threads=null num_io_threads=null num_bam_file_handles=null read_group_black_list=null pedigree=[] pedigreeString=[] pedigreeValidationType=STRICT allow_intervals_with_unindexed_bam=false generateShadowBCF=false logging_level=INFO log_to_file=null help=false genotype_likelihoods_model=SNP p_nonref_model=EXACT pcr_error_rate=1.0E-4 noSLOD=false annotateNDA=false min_base_quality_score=17 max_deletion_fraction=0.05 cap_max_alternate_alleles_for_indels=false min_indel_count_for_genotyping=5 min_indel_fraction_per_sample=0.25 indel_heterozygosity=1.25E-4 indelGapContinuationPenalty=10 indelGapOpenPenalty=45 indelHaplotypeSize=80 noBandedIndel=false indelDebug=false ignoreSNPAlleles=false allReadsSP=false ignoreLaneInfo=false reference_sample_calls=(RodBinding name= source=UNBOUND) reference_sample_name=null sample_ploidy=2 min_quality_score=1 max_quality_score=40 site_quality_prior=20 min_power_threshold_for_calling=0.95 min_reference_depth=100 exclude_filtered_reference_sites=false heterozygosity=0.001 genotyping_mode=DISCOVERY output_mode=EMIT_VARIANTS_ONLY standard_min_confidence_threshold_for_calling=30.0 standard_min_confidence_threshold_for_emitting=10.0 alleles=(RodBinding name= source=UNBOUND) max_alternate_alleles=3 dbsnp=(RodBinding name=dbsnp source=/srv/gs1/projects/cherry/giltae/hugeseq_data/snp.d/chrI_SNPs.vcf) comp=[] out=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub no_cmdline_in_header=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub sites_only=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub bcf=org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub debug_file=null metrics_file=null annotation=[DepthOfCoverage, MappingQualityZero] excludeAnnotation=[] filter_mismatching_base_and_quals=false"
# ##contig=<ID=chrI,length=230218>
# ##contig=<ID=chrII,length=813184>
# ##contig=<ID=chrIII,length=316620>
# ##contig=<ID=chrIV,length=1531933>
# ##contig=<ID=chrV,length=576874>
# ##contig=<ID=chrVI,length=270161>
# ##contig=<ID=chrVII,length=1090940>
# ##contig=<ID=chrVIII,length=562643>
# ##contig=<ID=chrIX,length=439888>
# ##contig=<ID=chrX,length=745751>
# ##contig=<ID=chrXI,length=666816>
# ##contig=<ID=chrXII,length=1078177>
# ##contig=<ID=chrXIII,length=924431>
# ##contig=<ID=chrXIV,length=784333>
# ##contig=<ID=chrXV,length=1091291>
# ##contig=<ID=chrXVI,length=948066>
# ##contig=<ID=chrM,length=85779>
# ##reference=file:///srv/gs1/projects/cherry/sacCer3/sacCer3.fa
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BY4741
# chrI	12529	.	A	G	794.01	PASS	AC=1;AF=0.500;AN=2;BaseQRankSum=-3.132;DP=155;Dels=0.00;FS=16.877;HaplotypeScore=3.4557;MLEAC=1;MLEAF=0.500;MQ=57.86;MQ0=0;MQRankSum=-0.583;QD=5.12;ReadPosRankSum=-3.218;SB=-9.501e+01	GT:AD:DP:GQ:PL	0/1:116,39:155:99:824,0,3691
# chrI	12620	.	T	C	1223.01	PASS	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.964;DP=321;Dels=0.00;FS=3.420;HaplotypeScore=1.6208;MLEAC=1;MLEAF=0.500;MQ=51.74;MQ0=3;MQRankSum=-1.011;QD=3.81;ReadPosRankSum=-1.500;SB=-8.360e+02	GT:AD:DP:GQ:PL	0/1:227,94:321:99:1253,0,6046

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
$CHR{"chrM"} = "chrM";

$line = "";
while($line = <IN>) {
	chomp($line);

	if($line =~ "#") {
		if($line =~ "##contig") {
			@array = split(/ID=/, $line);
			@subsplit = split(/,/, $array[1]);
			$subsplit[0] = $CHR{$subsplit[0]};
			print OUT $array[0] . "ID=" . $subsplit[0] . "," . $subsplit[1] . "\n";
		} else {
			print OUT $line . "\n";
		}
		next
	}
	@array = split(/\t/, $line);
	$array[0] = $CHR{$array[0]};
	# Print out gene info
	print OUT join("\t", @array),"\n";

}
close IN;
close OUT;
