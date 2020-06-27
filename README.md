# Bioinformatic Pipeline for Processing lacI Library

## Pipeline Overview

- This pipeline performs initial QC of the input PacBio sequencing data.
- The required input is a g-zipped fastq file contining PacBio CCS reads.
- Uses the lacI sequence to determine read direction and reverse complements reverse reads.
- Extracts target sequences, e.g. lacI and barcodes from reads.
- Generates QC stats for extracted reads.

## Pipeline Environment
- Snakemake is used to compose the bioinformatic pipeline, https://snakemake.readthedocs.io/en/stable/.  
- Snakemake and pipeline dependencies can be install using conda (https://www.anaconda.com/distribution/) and bioconda (https://bioconda.github.io/).
- The pipeline is defined in the `Snakefile` and adapter sequences are defined in the `targets.csv` config file. 

## Environment Setup
- Install conda with python 3.7 if not already installed on the system, https://docs.anaconda.com/anaconda/install/, miniconda as well as the full anaconda install will work. 
- Generate a conda environment for running the pipeline, `conda env create -n lacI_ccs --file environment.yaml`.


## Running pipeline
- Activate conda environment `conda activate lacI_ccs`.
- Run pipeline using `snakemake -j [#]`. The `-j` parameter takes the number of threads or parallel jobs snakemake executes simultaneously. (I like to also include `-p`, which prints the command snakemake executes at each step.) The dry run snakemake argument, `-n`, can be used to see which pipeline steps will be run.

## Modifying pipeline
- Pipeline input files and output directory are defined in the config file, 'config.yaml'. 
- The adapter sequences are defined and can be updated in the `targets.csv` file. 
- Snakemake works similar to make in that it only reruns parts of the analysis pipeline when an input file is changed or the output file is missing. After updating the `targets.csv` file you will want to remove any pipeline output files you want to replace.  

## Test data
- The file, '/data_0/raw/test.fastq.gz', contains a suitable test dataset. To run the pipeline with another dataset, change the output directory (outdir) and input file (infq) in the config.yaml file.