#!/usr/bin/env python3
import sys
import re

floatregex="(^.+) \(([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\)"

def main(argv):
	if (len(argv)<1):
		print("usage rnafold2csv.py <inputfile>")
		print("Converts rnafold result file into an .csv to STDOUT")
		exit(1)
	inputfile = argv[0]
	with open(inputfile) as ifile:
		print('"comment","sequence","secstructure","mfe"')
		for commentline in ifile:
			comment=commentline.replace("\n","").replace("\r","").replace(">","")
			sequence=next(ifile).replace("\n","").replace("\r","")
			secandmfeline=next(ifile).replace("\n","").replace("\r","")
			mferesult=re.search(floatregex,secandmfeline)
			if mferesult:
				mfe=mferesult.group(2)
				secstruct=mferesult.group(1)
			else:
				mfe=0
				secstruct=""
			print(f'"{comment}","{sequence}","{secstruct}",{mfe}')
if __name__ == "__main__":
   main(sys.argv[1:])
