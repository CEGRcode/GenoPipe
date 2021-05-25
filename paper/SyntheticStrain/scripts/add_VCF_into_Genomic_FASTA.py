import sys
import argparse
import random


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
	parser = argparse.ArgumentParser(description='Inclusive tool for SimulateStrains that helps to test Genopipe\'s StrainID tool')
	parser.add_argument('-f','--fasta', metavar='ref_genome_fn', required=True, help='a FASTA file of reference genome in FASTA format')
	parser.add_argument('-v','--vcf', metavar='strain_vcf_fn', required=True,help='the vcf file with which to make the genome ')
	parser.add_argument('-s','--size', default=60, type=int, help='a format value, the number of bp to print per line in output FASTA')
	args = parser.parse_args()
	return(args)

def load_vcf(vcf_fn):
	'''Load VCF file into a dictionary such that D[chr][pos] returns (ref,alt)'''
#-       CTG     0       PASS    .
#chr10   124598545       .       AGA     -       0       PASS    .
#chr19   37733623        .       TAAT    -       0       PASS    .
#chr19   49079290        .       C       -       0       PASS    .
#chr2    207173230       .       T       -       0       PASS    .
#chr5    23510054        .       C       -
	chr2variant = {}
	vcf_reader = open(vcf_fn,'r')
	for line in vcf_reader:
		if( line.find('#')==0 ):
			continue   # skip header lines
		tokens = line.strip().split()
		c = tokens[0]
		chr2variant.setdefault(c,[])
		pos = int(tokens[1])
		ref = tokens[3]
		alt = tokens[4]
		chr2variant[c].append( (pos,ref,alt) )
	vcf_reader.close()
	return( chr2variant )

def write_seq(seq,size):
	'''Write line width and return leftover seq'''
	while(size<len(seq)):
		sys.stdout.write('%s\n' % seq[:size] )
		seq = seq[size:]
	return( seq )


def add_variants(seq, pos2variant):
	'''Write sequence to standard out with variants included'''
	sys.stderr.write('len before: %i\n' % len(seq))

	seq = seq.upper()

	# vcf pos has 1st base at index 1

	# Update sequence with each variant
	for variant in reversed(sorted(pos2variant)):
		pos = variant[0]
		v_ref = variant[1]
		v_alt = variant[2]

		ref_tok = v_ref.split(',')
		alt_tok = v_alt.split(',')

		# Skip var if multiple refs specified
		if( len(ref_tok)>1 ):
			sys.stderr.write('Skip multi-reference entries: @%i : %s > %s\n' % (pos,v_ref,v_alt) )
#			# Check that all reference options are valid
#			if(all(nt in NT for nt in ref_tok)):
#				# Check reference string
#				fa_start = pos-1
#				fa_len = len(ref_tok)
#				if( seq[fa_start,fa_start+fa_len]!=v_ref ):
#					sys.stderr.write('VCF Reference does not match FASTA Reference: @pos %i, %s != %s\n' % (pos, seq[:pos-1], v_ref))
#					continue
#			# [ATCG]+ >> -
#			if(v_alt=='-'):
#				seq = seq[:fa_start]
			continue

		# Take first alt var if multiple available
		if( len(alt_tok)>1 ):
			sys.stderr.write('Multiple alternates available, pick first one @%i : %s > %s\n' % (pos,v_ref,v_alt) )
			v_alt = alt_tok[0]

		v_ref = v_ref.upper()
		v_alt = v_alt.upper()

		ref_len = len(v_ref)
		if(v_ref=='-'):
			# - >> [ATCG]+
			if( all(nt in NT for nt in list(v_alt)) ):
				seq = seq[:pos-1] + v_alt + seq[pos-1+ref_len:]
			# - >> unaccounted for allele
			else:
				sys.stderr.write('(CASE UNACCOUNTED) Skip variant @%i : %s > %s\n' % (pos,v_ref,v_alt) )
		elif( all(nt in NT for nt in list(v_ref)) ):
			# Check reference string
			if( seq[pos-1:pos-1+ref_len]!=v_ref ):
				sys.stderr.write('VCF Reference does not match FASTA Reference: @pos %i, %s != %s\n' % (pos, seq[pos-1:pos-1+ref_len], v_ref))
				continue
			# [ATCG]+ >> -
			if( v_alt=='-' ):
				seq = seq[:pos-1] + seq[pos-1+ref_len:]
			# [ATCG]+ >> [ATCG]+
			elif( all(nt in NT for nt in list(v_alt)) ):
				seq = seq[:pos-1] + v_alt + seq[pos-1+ref_len:]
			# [ATCG]+ >> unaccounted for allele
			else:
				sys.stderr.write('(CASE UNACCOUNTED) Skip variant @%i : %s > %s\n' % (pos,v_ref,v_alt) )
		else:
			sys.stderr.write('(CASE UNACCOUNTED) Skip variant @%i : %s > %s\n' % (pos,v_ref,v_alt) )
			#fill in code for handling more complex vcf entries
			# for now, (1) skipping and using ref
			# add	   (2) randomly picking one (when 2 alts present)
			# add	   (3) encode with degenerate codings for translation during random seq generation

	sys.stderr.write('len after: %i\n' % len(seq))
	return seq

def format_print(seq, size):
	'''Write sequence to STDOUT updated to write once'''
	seq_w_newlines = '\n'.join(seq[i:i+size] for i in range(0, len(seq), size))
	sys.stdout.write( seq_w_newlines+'\n' )

def format_print1(seq, size):
	'''Write sequence to STDOUT by the designated character count line size'''
	while( size<len(seq) ):
		sys.stdout.write('%s\n' % seq[:size] )
		seq = seq[size:]
	sys.stdout.write('%s\n' % seq)

# Main program which takes in input parameters
if __name__ == '__main__':
	'''Given a ref sequence and a corresponding VCF file, generate the strain genome.
	Note that the accounting for multiple variants is more robust to simpler variants'''
	args = getParams()

	ref_reader = open(args.fasta, 'r')
	sys.stderr.write('Reference genome: %s\n' % args.fasta)

	chr2vlist = load_vcf(args.vcf)
	sys.stderr.write('Loaded variants: %s\n' % args.vcf)
	seq = ''
	variants = None
	raw_line = None
	# For each line in genome FASTA
	while(raw_line!=""):
		# Parse next line
		raw_line = ref_reader.readline().strip()
		# Handle header lines:
		if(raw_line.find('>')==0):
			# Write new seq of current chr before parsing new chr info
			if(variants != None):
				# Add variants for current chr seq
				new_seq = add_variants(seq,variants)
				# Write this new sequence formatted by standard seq size
				format_print(new_seq, args.size )
			# Write new seq FASTA header
			sys.stdout.write(raw_line + "\n")
			# Reset variables
			seq = ''
			chr_key = raw_line[1:]
			variants = chr2vlist.get(chr_key,[])
			variants.sort()
			sys.stderr.write('%s: (%i variants)\n' % (chr_key,len(variants)) )
			continue
		seq += raw_line
	new_seq = add_variants(seq,variants)
	format_print(new_seq, args.size )
