# This is the main file for Snakemake, it defines the workflow steps, which need to be executed.https://snakemake.readthedocs.io/en/v5.3.0/index.html


# The configuration file to load
configfile: "config.yaml"

noofsplits=2
splitindices=['%07d'%i for i in range(0,noofsplits)];


# Paths
basedir = config["datadir"]
runshuffleinstall = basedir+"/installedshuffle"
seqshuffled = basedir+"/shuffled.fst"

shufflelist=[20,100,200,1000]
shufflemethods=["m","d","z","f"]
# m, mononucleotide shuffling; d, dinucleotide shuffling; z, zero-order markov model; f, first-order markov model

wildcard_constraints:
    index="\d+"

inputgroups=["real_izmar","pseudo_izmar"]
# Real miRNA has to contain "real"


localrules: presentation, arff, splitfasta, mergecsv, mergefinalcsv, fasta2csv, joincsv


#
# Run rnafold
#
rule fold:
	input:
		basedir+"/{inputgroup}/split/{inputgroup}.fasta_chunk_{index}"
	output:
		basedir+"/{inputgroup}/fold/{index}.fold"
	conda: 
		"envs/rnafold.yaml"
	shadow:
		"shallow"
	group:
		"fold"
	shell:
		"RNAfold --noPS -p < {input} > {output}"

#
# Run stanley genRNAStats.pl
#
rule stanleyRNAstats:
	input:
		basedir+"/{inputgroup}/split/{inputgroup}.fasta_chunk_{index}"
	output:
		stats=basedir+"/{inputgroup}/stanley/{index}.stats"
	group:
		"stanleyfeatures"
	shell:
		"perl scripts/shuffle/genRNAStats.pl -i {input} -o {output.stats}"
#
# Parse the stanley features
#
rule parsestnlyfeatures:
	input:
		shuffledstatfiles=rules.stanleyRNAstats.output.stats
	output:
		basedir+"/{inputgroup}/datasplit/{index}.stnley.csv"
	group:
		"stanleyfeatures"
	shell:
		"cp {input} {output}"
#
# Shuffle the sequences
#
rule snuffleshuffel:
	input:
		basedir+"/{inputgroup}/split/{inputgroup}.fasta_chunk_{index}"
	output:
		basedir+"/{inputgroup}/stanley/shuffled/{index}-{method}-{shuffles}.fasta"
	params:
		method="{method}",
		shuffles="{shuffles}"
	group:
		"shuffle"
	shell:
		"perl scripts/shuffle/genRandomRNA.pl -n {params.shuffles} -m {params.method} < {input} > {output}"
#
# Fold the shuffled sequences
#
rule foldshuffled:
	input:
		rules.snuffleshuffel.output
	output:
		basedir+"/{inputgroup}/stanley/shuffled/{index}-{method}-{shuffles}.fold"
	group:
		"shuffle"
	shell:
		"RNAfold --noPS < {input} > {output}"
#
# Compute Stanleys features from the shuffled sequences
#	
rule stnlyRandfeatures:
	input:
		unshuffled=rules.fold.output,
		shuffled=rules.foldshuffled.output,
	output:
		stats=basedir+"/{inputgroup}/stanley/{index}-{method}-{shuffles}.stats"
	params:
		shuffles="{shuffles}"
	group:
		"shuffle"
	shell:
		"perl scripts/shuffle/genRNARandomStats.pl -n {params.shuffles} -i {input.shuffled} -m {input.unshuffled} -o {output.stats}"

#
# Parse the stanley randfold features
#
rule parsestnlyRandfeatures:
	input:
		shuffledstatfiles=rules.stnlyRandfeatures.output.stats
	output:
		stats=basedir+"/{inputgroup}/stanley/{index}-{method}-{shuffles}.csv"
	group:
		"shuffle"
	shell:
		"cp {input} {output}"

#
# Run RNAspectral
#
rule RNAspectral:
	input:
		rules.fold.output
	output:
		basedir+"/{inputgroup}/stanley/{index}.spectral"
	group:
		"spectral"
	shell:
		# Clean the additional information
		"grep --invert-match '[]}}]$\| frequ' {input} | scripts/shuffle/RNAspectral.exe > {output}"
#
# Parse RNAspectral output
#
rule parseRNAspectral:
	input:
		shuffledstatfiles=rules.RNAspectral.output
	output:
		basedir+"/{inputgroup}/datasplit/{index}.spectral.csv"
	group:
		"spectral"
	shell:
		"cp {input} {output}"

#
# Parse rnafold
#
rule parsernafold:
	input:
		rules.fold.output
	output:
		basedir+"/{inputgroup}/datasplit/{index}.fold.csv"
	group:
		"fold"
	script:
		"scripts/rnafold2csv/rnafold2csv.py"
#
# Split fasta file
#
rule splitfasta:
	input:
		basedir+"/{inputgroup}.fasta"
	output:
		expand(basedir+"/{{inputgroup}}/split/{{inputgroup}}.fasta_chunk_{index}",index=splitindices)
	params:
		splits=noofsplits,
		outputdir=directory(basedir+"/{inputgroup}/split/")
	shell:
		"fastasplit -f {input} -o {params.outputdir} -c {params.splits}"

#
# Convert the given fasta chunk into a csv file
#
rule fasta2csv:
	input:
		basedir+"/{inputgroup}/split/{inputgroup}.fasta_chunk_{index}"
	output:
		basedir+"/{inputgroup}/datasplit/{index}.seq.csv"
	params:
		realmarker="real"
	script:
		"scripts/fasta2csv/fasta2csv.R"
#
# Calculates features from the fold csv file
#
rule derviedcsv:
	input:
		rules.parsernafold.output
	output:
		basedir+"/{inputgroup}/datasplit/{index}.derived.csv"
	group:
		"fold"
	script:
		"scripts/features_derived/features_derived.R"
#
# Join the calculated .csv files
#
rule joincsv:
	input:
		expand(basedir+"/{{inputgroup}}/datasplit/{{index}}.{type}.csv",type=["fold","seq","derived","stnley","spectral"]),
		expand(basedir+"/{{inputgroup}}/stanley/{{index}}-{method}-{shuffles}.csv",method=shufflemethods,shuffles=shufflelist) 
	output:
		basedir+"/{inputgroup}/split-{index}.csv"
	script:
		"scripts/csvmerge/csvmerge.R"

#
# Merge the .csv files from the sets
#
rule mergecsv:
	input:
		csvs=expand(basedir+"/{{inputgroup}}/split-{index}.csv",index=splitindices)
	output:
		csv=basedir+"/{inputgroup}/combined.csv"
	script:
		"scripts/concatenateCsvs/concatenateCsvs.R"
#
# Merge the generated .csv files
#
rule mergefinalcsv:
        input:
                csvs=expand(rules.mergecsv.output.csv,inputgroup=inputgroups)
        output:
                csv=basedir+"/all.csv"
        script:
                "scripts/concatenateCsvs/concatenateCsvs.R"

#
# Generate .arff for Weka
#
rule arff:
	input:
		rules.mergefinalcsv.output.csv
	output:
		basedir+"/all.arff"
	script:
		"scripts/csv2arff/csv2arff.R"


#
# Generate the project presentation
#
rule presentation:
	input:
		template="presentation/template.pptx",
		finalcsv=rules.mergefinalcsv.output.csv
	output:
		presentation=basedir+"/presentation.pptx"
	script:
		"presentation/projectpresentation.Rmd"



#rule dustmasker:
#    input: realhairpins
 #   output: realhairpins+"-dust"
  #  shell: "dustmasker -in {output} -outfmt fasta -out {output} -level 15"

#rule installPerlShuffle:
#	output: runshuffleinstall
#	shell: "PERL_MM_USE_DEFAULT=1 cpan Algorithm::Numerical::Shuffle > {output}"
#

