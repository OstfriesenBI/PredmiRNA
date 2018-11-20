# This is the main file for Snakemake, it defines the workflow steps, which need to be executed.https://snakemake.readthedocs.io/en/v5.3.0/index.html


# The configuration file to load
configfile: "config.yaml"

basedir = config["datadir"]
SAMPLES = ["1", "2"]
realhairpins = basedir+"/mirbasehairpin.fst"


# Just to test it for now:

# SnakeMake creates directories automagically, if they are defined in the input/outputs
rule all:
    input: realhairpins



rule downloadhairpinlong:
    output: realhairpins+"-long"
    shell: "curl ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz | gunzip > {output}"

rule shorthairpin:
    input: realhairpins+"-long"
    output: realhairpins
    shell: "head -n101 {input} >{output}"
