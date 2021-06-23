import sys
import argparse
import random
import pysam

degenerate_codes = {
	'C,G':'S','G,C':'S',
	'A,C':'M','C,A':'M',
	'A,G':'R','G,A':'R',
	'A,T':'W','T,A':'W',
	'C,T':'Y','T,C':'Y',
	'G,T':'K','T,G':'K' }
NT = ['A','T','C','G']
sgd2chr = { 'chrI':'chr1',
	'chrII':'chr2', 'chrIII':'chr3',
	'chrIV':'chr4','chrV':'chr5',
	'chrVI':'chr6','chrVII':'chr7',
	'chrVIII':'chr8','chrIX':'chr9',
	'chrX':'chr10','chrXI':'chr11',
	'chrXII':'chr12','chrXIII':'chr13',
	'chrXIV':'chr14','chrXV':'chr15',
	'chrXVI':'chr16','chrmt':'chrM' }

# since VCF files seem to be pretty small, we can probably load it all into a dictionary,
# then randomly sample

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Script to incorportate VCF alternate alleles into a reference genomic FASTA.')
	parser.add_argument('-f','--fasta', metavar='ref_genome_fn', required=True, help='a FASTA file of reference genome in FASTA format')
	parser.add_argument('-v','--vcf', metavar='strain_vcf_fn', required=True,help='the vcf file with which to make the genome ')
	parser.add_argument('-s','--size', default=60, type=int, help='a format value, the number of bp to print per line in output FASTA')
	args = parser.parse_args()
	return(args)

def add_variants(seq, var):
	'''Write sequence to standard out with variants included. Expect input list of VariantRecord objects to be reverse sorted by reference start coordinate. Assumes no overlap between start and stop intervals of any two variants in VCF.'''
	START_LEN = len(seq)
	seq = seq.upper()
	# Update sequence with each variant
	for VCFrecord in var:
		# Get first alt allele length
		alen = len(VCFrecord.alts[0])
		# Skip record if multiple alt alleles specified
		if(len(VCFrecord.alts)>1):
			sys.stderr.write("Skip multi-alternate allele entries: @%i : %s > %s\n" % (VCFrecord.pos,VCFrecord.ref,"|".join(VCFrecord.alts)))
			continue
		## Substitution handled
		elif(VCFrecord.ref.upper() in NT and VCFrecord.alts[0].upper() in NT):
			# Check expected ref
			if(seq[VCFrecord.start]!=VCFrecord.ref):
				sys.stderr.write("(!) Skip reference mismatch: expected %s != found %s\n" % (VCFrecord.ref,seq[VCFrecord.start]))
				continue
			sys.stderr.write("(*) Process Substitution: %s\n" % (str(VCFrecord).strip()))
			sys.stderr.write("Length before substitution:\t%i\t" % len(seq))
			seq = seq[:VCFrecord.start] + VCFrecord.alts[0] + seq[VCFrecord.start+1:]
			sys.stderr.write("...and after:\t%i\n" % len(seq))
			continue
		## Insertion handled
		elif(VCFrecord.ref=="-"):
			if(alen>1):
				sys.stderr.write("(!) Skip Insertions >1nt long: %s\n" % (str(VCFrecord).strip()))
				continue
			sys.stderr.write("(*) Process Insertion: %s\n" % (str(VCFrecord).strip()))
			sys.stderr.write("Length before inserton:\t%i\t" % len(seq))
			seq = seq[:VCFrecord.start] + VCFrecord.alts[0] + seq[VCFrecord.start:]
			sys.stderr.write("...and after:\t%i\n" % len(seq))
			continue
		## Deletion handled
		elif(VCFrecord.alts[0]=="-" or VCFrecord.rlen>alen):
			if(VCFrecord.rlen>1):
				sys.stderr.write("(!) Skip Deletions >1nt long: %s\n" % (str(VCFrecord).strip()))
				continue
			# Check expected ref
			if(seq[VCFrecord.start:VCFrecord.stop]!=VCFrecord.ref):
				sys.stderr.write("(!) Skip reference mismatch: expected %s != found %s\n" % (VCFrecord.ref,seq[VCFrecord.start:VCFrecord.start+VCFrecord.rlen]))
				continue
			sys.stderr.write("(*) Process Deletion: %s\n" % (str(VCFrecord).strip()))
			sys.stderr.write("Length before deletion:\t%i\t" % len(seq))
			# seq = seq[:VCFrecord.start] + VCFrecord.alts[0] + seq[VCFrecord.start+VCFrecord.rlen:]
			seq = seq[:VCFrecord.start] + seq[VCFrecord.start+VCFrecord.rlen:]
			sys.stderr.write("...and after:\t%i\n" % len(seq))
			continue
		## Invalid VCF entry handled
		else:
			sys.stderr.write("(!) INVALID OR UNHANDLED VCF ENTRY: %s\n" % str(VCFrecord).strip())
			continue


	sys.stderr.write('Total length: %i > %i\n' % (START_LEN, len(seq)))
	return seq

def format_print(seq, size):
	'''Write sequence to STDOUT updated to write once'''
	seq_w_newlines = '\n'.join(seq[i:i+size] for i in range(0, len(seq), size))
	sys.stdout.write( seq_w_newlines+'\n' )


# Main program which takes in input parameters
if __name__ == '__main__':
	'''Given a ref sequence and a corresponding VCF file, generate the strain genome.
	Note that the accounting for multiple variants is more robust to simpler variants'''
	args = getParams()

	sys.stderr.write('Reference genome: %s\n' % args.fasta)
	FASTA = pysam.FastaFile(args.fasta)

	# chr2vlist = load_vcf(args.vcf)
	sys.stderr.write('Variants File: %s\n' % args.vcf)
	vcffile = pysam.VariantFile(args.vcf)
	all_variants = list(vcffile.fetch())

	for CHR in FASTA.references:
		sys.stderr.write(CHR)
		variants = sorted([v for v in all_variants if v.contig==CHR], key = lambda x: x.start, reverse=True)
		new_seq = add_variants(FASTA.fetch(CHR),variants)
		sys.stdout.write(">%s\n" % (CHR))
		format_print(new_seq, args.size)
