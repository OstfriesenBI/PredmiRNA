#!/bin/bash
if [ $# -ne 4 ]; then
    echo "Usage: rnafold.sh <inputfile> <outputfile> <alg flag: -s/-d/-m> <noofshuffles>"
    exit 1
fi
input=$1
output=$2
flag=$3
shuffles=$4
tempfile=$(mktemp)

echo "\"comment\",\"randfoldp\"" > $2
sed "s/ /;/g" $1 > $tempfile
randfold $flag $tempfile $shuffles | sed "s/;/ /g" | awk -F"\t" '{print "\"" $1 "\",\"" $3 "\""}' >> $2

rm $tempfile
