import sys
import re
import argparse

def getParams():
	'''Parse parameters from the command line'''
	parser = argparse.ArgumentParser(description='Parse runtimes from system and convert to seconds for plotting purposes.')
	parser.add_argument('-i','--input-times', metavar='grep_time', required=True, help='the file "grepping" out the runtimes, `grep "real" input_dir/*.time`')
	args = parser.parse_args()
	return(args)

def parse_time(time_string):
	parsed_time = re.findall("([0-9]+)m([\.0-9]+)s", time_string)
	seconds = float(parsed_time[0][1])
	minutes = int(parsed_time[0][0])
	sys.stdout.write("%f\n" % (seconds + 60.0*minutes))

if __name__ == "__main__":
	args = getParams()
	reader = open(args.input_times)
	for line in reader:
		tokens = line.strip().split("\t")
		parse_time(tokens[1])
	reader.close()
