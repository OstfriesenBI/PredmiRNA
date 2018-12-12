# This is the main file for Snakemake, it defines the workflow steps, which need to be executed.https://snakemake.readthedocs.io/en/v5.3.0/index.html


# The configuration file to load
configfile: "config.yaml"

basedir = config["datadir"]
SAMPLES = ["1", "2"]
realhairpins = basedir+"/mirbasehairpin.fst"
rnafoldout = basedir+"/folded.txt"
rnafoldps = basedir+"/rnaps/"
rnafoldcsv = basedir+"/folded.csv"
# Just to test it for now:

# SnakeMake creates directories automagically, if they are defined in the input/outputs
rule all:
    input: rnafoldcsv


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
    shell: "./scripts/rnafold2csv/rnafold2csv.py {input} > {output}"
