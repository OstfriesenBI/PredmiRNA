#!/usr/bin/env python3

import csv
import re

filetoread="data/pseudo_izmir/datasplit/0000001.fold.csv"

with open(filetoread,'r') as inputcsv:
	seqs=[[row["sequence_nodm"],row["mfesecstructure"]] for row in csv.DictReader(inputcsv)]

p=re.compile("\\.{4,}")

#loops=[ [m.start(), m.end()] for m in p.finditer(seqi[1]) for seqi in seqs]
loops=[[ i for s in p.finditer(seqi[1]) for i in range(s.start(),s.end())] for seqi in seqs]

print(seqs[0][1])
print(loops[0])
