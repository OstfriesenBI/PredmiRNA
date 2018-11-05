# This is the main file for Snakemake, it defines the workflow steps, which need to be executed.https://snakemake.readthedocs.io/en/v5.3.0/index.html


# The configuration file to load
configfile: "config.yaml"

basedir = config["datadir"]
SAMPLES = ["1", "2"]


# Just to test it for now:

# SnakeMake creates directories automagically, if they are defined in the input/outputs
rule all:
    input: expand("data/input{sample}.txt",sample=SAMPLES)

rule compute1:
    output: "data/input1.txt"
    shell: "touch data/input1.txt"

rule compute2:
    conda: "envs/rnafold.yaml"
    output: "data/input2.txt"
    shell: "RNAfold -h > {output}"
