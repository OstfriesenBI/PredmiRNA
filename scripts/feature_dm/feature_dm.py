#!/usr/bin/env python3
import sys
import re

def cleanstring(toclean):
	return toclean.replace("\n","").replace("\r","")


def n_lower_chars(string):
    return sum(1 for c in string if c.islower())

def main(argv):
	if (len(argv)<1):
		print("usage feature_dm.py <inputfile>",file=sys.stderr)
		print("Reads a FASTA file output by dustmasker and wirtes a csv file with the comment and the low complexity percentage.",file=sys.stderr)
		exit(1)
	inputfile = argv[0]
	with open(inputfile) as ifile:
		print('"comment","dm"')
		comment = ''
		sequence = ''
		for line in ifile:
			if(line.startswith(">")):
				if(comment != ''):
					dm=n_lower_chars(sequence)/len(sequence)*100
					print(f'"{comment}","{dm}"')
				comment = cleanstring(line).replace(">","") # Just remove the first char, this may go badly, you watch...
				sequence = ''
			else:
				sequence = sequence + cleanstring(line)
		print(f'"{comment}","{dm}')
if __name__ == "__main__":
   main(sys.argv[1:])
