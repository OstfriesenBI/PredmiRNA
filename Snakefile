# This is the main file for Snakemake, it defines the workflow steps, which need to be executed.https://snakemake.readthedocs.io/en/v5.3.0/index.html


# The configuration file to load
configfile: "config.yaml"

noofsplits=1
splitindices=['%07d'%i for i in range(0,noofsplits)];


# Paths
basedir = config["datadir"]
realhairpins = basedir+"/mirbasehairpin.fst"
rnafoldout = basedir+"/folded.txt"
rnafoldps = basedir+"/rnaps/"
rnafoldcsv = basedir+"/folded.csv"
randfoldcsv = basedir+"/randfold.csv"
runshuffleinstall = basedir+"/installedshuffle"
seqshuffled = basedir+"/shuffled.fst"

wildcard_constraints:
    index="\d+"

inputgroups=["real_izmar","pseudo_izmar"]
# Real miRNA has to contain "real"

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
	shell:
		"RNAfold --noPS -p < {input} > {output}"
#
# Parse rnafold
#
rule parsernafold:
	input:
		rules.fold.output
	output:
		basedir+"/{inputgroup}/datasplit/{index}.fold.csv"
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
# Join the calculated .csv files
#
rule joincsv:
	input:
		expand(basedir+"/{{inputgroup}}/datasplit/{{index}}.{type}.csv",type=["fold","seq"]) 
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


rule rnafold2csv:
    input: rnafoldout
    output: rnafoldcsv
    shell: "scripts/rnafold2csv/rnafold2csv.py {input} > {output}"

rule dustmasker:
    input: realhairpins
    output: realhairpins+"-dust"
    shell: "dustmasker -in {output} -outfmt fasta -out {output} -level 15"

rule installPerlShuffle:
	output: runshuffleinstall
	shell: "PERL_MM_USE_DEFAULT=1 cpan Algorithm::Numerical::Shuffle > {output}"

rule shuffleSeq:
	input:
		seq=realhairpins
	output: 
		seqshuffled
	shell: "perl scripts/shuffle/genRandomRNA.pl -n 200 -m m < {input.seq} > {output}"
      # m d z oder f


#For RNAspectral: grep --invert-match '[]}]$\| frequ' folded.txt | ../scripts/shuffle/RNAspectral.exe
