# This is the main file for Snakemake, it defines the workflow steps, which need to be executed.https://snakemake.readthedocs.io/en/v5.3.0/index.html


# The configuration file to load
configfile: "config.yaml"

splitindices=[i for i in range(0,100)];


# Paths
basedir = config["datadir"]
realhairpins = basedir+"/mirbasehairpin.fst"
rnafoldout = basedir+"/folded.txt"
rnafoldps = basedir+"/rnaps/"
rnafoldcsv = basedir+"/folded.csv"
randfoldcsv = basedir+"/randfold.csv"
runshuffleinstall = basedir+"/installedshuffle"
seqshuffled = basedir+"/shuffled.fst"


inputgroups=["real_izmar","pseudo_izmar"]


#
# Join the calculated .csv files
#
rule joincsv:
	output:
		basedir+"/{inputgroup}/split-{index}.csv"
	shell:
		"touch {output}"

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


rule downloadhairpinlong:
    output: realhairpins+"-long"
    shell: "curl ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz | gunzip > {output}"

rule shorthairpin:
    input: realhairpins+"-long"
    output: realhairpins
    shell: "head -n101 {input} >{output}"

rule rnafold:
    output: rnafoldout
    input: realhairpins
    conda: "envs/rnafold.yaml"
    shell: "RNAfold --noPS -p < {input} > {output}"

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
