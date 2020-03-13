# LacI PacBio HiFi Analysis Pipeline

## Objective
Develop snakemake pipeline to process PacBio HiFi data for the LacI landscape manuscript. 

## Approach
- use seqkit and cutadapt for sequence manipulation, https://bioinf.shenwei.me/seqkit

## Identify forward and reverse seqs
- identify rev comp sequences: 
    - seqkit fish -f data/ref/lacI.fa data/raw/XTACK.ccs.head.fastq 2> data/processed/pipe_dev/lacI_fish.tsv
- generate id list for forward seqs
    - awk '{ if ($7== "+") {print $1}}' data/processed/pipe_dev/lacI_fish.tsv > pos_seqs.lst 
- generate id list for reverse seqs
    - awk '{ if ($7== "-") {print $1}}' data/processed/pipe_dev/lacI_fish.tsv > rev_seqs.lst 

### Split by strand
- forward seq fastq 
    - seqkit grep --pattern-file pos_seqs.lst {fastq} > forward.fastq
split by strand
    - seqkit grep --pattern-file rev_seqs.lst {fastq} > reverse.fastq

### Generate combined forward and reverse seq set
- reverse complement reverse sequences
    - seqkit seq -pv reverse.fastq > rev_comp.fastq
- combine forward and reverse complemented fastq
    - cat forward.fastq reverse.fastq > all_forward.fastq

## Filter seqs
- Length filter to avoid overly long and duplicate match issues when generating tables
- seqkit seq -m 5000 -M 6000 ????

### extract target sequences
- cutadapt for target sequences, use appropriate length filters

### QA fastqs
- seqkit fx2tab -n -l -g -i {fastq}

### get AA sequences
- seqkit translate 

### Make table
use fq2tab and join
https://bioinf.shenwei.me/csvtk/
conda install -c bioconda csvtk


## Steps to add for QA/QC
- Mapping to reference and stats (useful for error rates)