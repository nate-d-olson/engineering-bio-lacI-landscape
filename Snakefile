## Import dependencies
import pandas as pd
from snakemake.utils import min_version

### set minimum snakemake version
min_version("5.5.0")

## Parse target metadata file 
targets = pd.read_csv("targets.tsv").set_index("target", drop = False)

## Wildcard constraints
SEQGROUP=["forward", "reverse", "combined"]
wildcard_constraints:
    seq_group="|".join(SEQGROUP)


## Defining output directory 
outdir="data/processed/XTACK_20200108"

## Input fastq 
infq = "data/raw/XTACK_20200108_S64049_PL100138141A-1_A01.ccs.fastq.gz"

## All
rule all:
    input:
        outdir + "/XTACK_20200108.tsv" , 
        expand(outdir + "/{seq_group}_stats.tsv", 
            seq_group = SEQGROUP),
        expand(outdir + "/targets/{target}_stats.tsv",
            target = list(set(targets["target"]))),
        expand(outdir + "/targets/XTACK_2020010_{target}.tsv.gz",
            target = list(set(targets["target"])))

## Initial QC
rule ccs_qa:
    input: infq
    output: outdir + "/XTACK_20200108.tsv" 
    shell: "seqkit fx2tab -nlgiH {input} > {output}"

## Identify forward and reverse seqs
rule find_rev:
    input:
        ref="data/ref/lacI.fa",
        fq=infq
    output: outdir + "/lacI_fish.tsv"
    shell: "seqkit fish -j 14 -f {input.ref} {input.fq} 2> {output}"

## Extract forward seqs
rule forward_seqs:
    input: 
        fish_tbl=outdir + "/lacI_fish.tsv",
        fq=infq
    output: temp(outdir + "/forward_seqs.fastq")
    shell: """
        awk '{{ if ($7== "+") {{print $1}}}}' {input.fish_tbl} \
            | seqtk subseq {input.fq} - > {output}
    """

## Extract and reverse complement reverse seqs
rule reverse_seqs:
    input: 
        fish_tbl=outdir + "/lacI_fish.tsv",
        fq=infq
    output: temp(outdir + "/reverse_seqs.fastq")
    shell: """
        awk '{{ if ($7== "-") {{print $1}}}}' {input.fish_tbl} \
            | seqtk subseq {input.fq} - |
            seqtk seq -r - > {output}
    """

## Generate combined forward and reverse seq set
rule combine_seqs:
    input: 
        forward=outdir + "/forward_seqs.fastq", 
        reverse=outdir + "/reverse_seqs.fastq"
    output: outdir + "/combined_seqs.fastq"
    shell: "cat {input.forward} {input.reverse} > {output}"

rule qc_seqs:
    input: outdir + "/{seq_group}_seqs.fastq"
    output: outdir + "/{seq_group}_stats.tsv"
    shell: "seqkit fx2tab -nlgiH {input} > {output}"
## TODO - add set to extract reads without lacI

## Filter seqs
## - Length filter to avoid overly long and duplicate match issues when generating tables
## - seqkit seq -m 5000 -M 6000 ????

### extract target sequences
def get_adp_param(wildcards):
    adp_flag=targets.loc[wildcards.target, "adp_flag"]
    adapter_seq=targets.loc[wildcards.target, "adapter"]
    
    min_len=targets.loc[wildcards.target, "min_len"]
    max_len=targets.loc[wildcards.target, "max_len"]

    return(f"{adp_flag} {adapter_seq} -m {min_len} -M {max_len}")

rule extract_target:
    input: outdir + "/combined_seqs.fastq"
    output: 
        fq=outdir + "/targets/{target}.fastq.gz",
        log=outdir + "/targets/{target}.log"
    params: 
        adp_param=get_adp_param
    shell: """
        cutadapt --discard-untrimmed \
            {params.adp_param} \
            -o {output.fq} \
            {input} \
           > {output.log}
    """

rule qa_targets:
    input: outdir + "/targets/{target}.fastq.gz"
    output: outdir + "/targets/{target}_stats.tsv"
    shell: "seqkit fx2tab -nlgiH {input} > {output}"

### get AA sequences
# - seqkit translate 

### Make table
rule make_seq_tbls:
    input: outdir + "/targets/{target}.fastq.gz"
    output: outdir + "/targets/XTACK_2020010_{target}.tsv.gz"
    shell: "seqkit fx2tab -H {input} | cut -f 1,2 | gzip -c > {output}" 
