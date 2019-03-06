#!/usr/bin/env python3
import sys
import re

'''
Output element:
>comment no 1
UACACAUAUUACCACCGGUGAACUAUGCAAUUUUCUACCUUACCGGAGACAGAACUCUUCGA
.....(.((..(((((((((((((.((((((((((((........(((((...))))))))))))))))).)))))))))))))..)).)......... (-42.50)
.....(.((..(((((((((((((.((((((((((((({,,.....,,{{...}}})))))))))))))).)))))))))))))..)).}......... [-44.14]
.....(.((..(((((((((((((.(((((((((((((..........((...))..))))))))))))).)))))))))))))..)).)......... {-38.20 d=6.37}
 frequency of mfe structure in ensemble 0.0698212; ensemble diversity 8.65  
'''

floatregex=r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
regexmfeline=r'(^.+) \( *('+floatregex+')\)'#Group 1: MFE sec. structure, Group 2: MFE (in kcal/mol)
regexprobline=r'(^.+) \[ *('+floatregex+')\]' #Group 1: base pair probabilities, Group 2: EFE (in kcal/mol)
regexcentroidline=r'(^.+) { *('+floatregex+') d=('+floatregex+')}' #Group 1: centroid sec. structure, Group 2: centroid MFE (in kcal/mol), Group 3: distance to ensemble
regexfreqdiv=r'frequency of mfe structure in ensemble ('+floatregex+'); ensemble diversity ('+floatregex+')' # Group 1: frequency, Group 2: diversity
def cleanstring(toclean):
	return toclean.replace("\n","").replace("\r","")

def main(argv):
	if (len(argv)<2):
		print("usage rnafold2csv.py <inputfile> <outputfile>",file=sys.stderr)
		print("Converts rnafold result file into an .csv to STDOUT",file=sys.stderr)
		exit(1)
	inputfile = argv[0]
	outputfile = argv[1]
	with open(outputfile,'w') as ofile:
		with open(inputfile) as ifile:
			print('"comment","sequence","mfesecstructure","mfe","basepairprobs","efe","centroidsecstructure","centroidmfe","ensembledistance","freqmfestruct","ensemblediversity"',file=ofile)
			for commentline in ifile:
				# Read the required lines:
				comment=cleanstring(commentline).replace(">","") # Just remove the first char, this may go badly, you watch...
				sequence=cleanstring(next(ifile))
				secandmfeline=cleanstring(next(ifile))
				pairprobsefeline=cleanstring(next(ifile))
				centroidmfeensembledist=cleanstring(next(ifile))
				frequencyandensemblediv=cleanstring(next(ifile))
				# Run regex search
				# MFE
				mferesult=re.search(regexmfeline,secandmfeline)
				if mferesult:
					mfe=mferesult.group(2)
					msecstruct=mferesult.group(1)
				else:
					print("Failed reading mfe/secondary structure in this line: "+secandmfeline,file=sys.stderr)
					print(regexmfeline)
					exit(2)
				# EFE
				ensembleresult=re.search(regexprobline,pairprobsefeline)
				if ensembleresult:
					efe=ensembleresult.group(2)
					bpprob=ensembleresult.group(1)
				else:
					print("Failed reading efe/base pair probabilities in this line: "+pairprobsefeline,file=sys.stderr)
					exit(2)
				# Centroid
				centroidresult=re.search(regexcentroidline,centroidmfeensembledist)
				if centroidresult:
					ensembledistance=centroidresult.group(3)
					centroidmfe=centroidresult.group(2)
					centroidsecstructure=centroidresult.group(1)
				else:
					print("Failed reading the centroid structure in this line: "+centroidmfeensembledist,file=sys.stderr)
					exit(2)
				#Frequency, Diversity
				freqresult=re.search(regexfreqdiv,frequencyandensemblediv)
				if freqresult:
					diversity=freqresult.group(2)
					frequency=freqresult.group(1)
				else:
					print("Failed reading the frequency/diversity in this line: "+frequencyandensemblediv,file=sys.stderr)
					exit(2)
				print(f'"{comment}","{sequence}","{msecstruct}",{mfe},"{bpprob}",{efe},"{centroidsecstructure}",{centroidmfe},{ensembledistance},{frequency},{diversity}',file=ofile)
try:
    snakemake
except NameError:
    snakemake = None

if snakemake is not None:
	main([snakemake.input[0],snakemake.output[0]])
elif __name__ == "__main__":
   main(sys.argv[1:])
