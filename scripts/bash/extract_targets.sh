#!/usr/bin/bash
FQ=$1

## Parameters
## --rc check for reverse compliment matches
## -m minimum length
## -M maximum length after trimming
extract_seqs () {
    cutadapt \
        --rc --discard-untrimmed -j 14 \
        -m 1 -M 5500 \
        -a  ${adapter} \
        --info-file ${outroot}.info \
        -o ${outroot}.fastq $FQ \
        > ${outroot}.log
}

# ## Barcode 1 ###################################################################
# adapter=ATCGGTGAGCCCGGGCTGTCGGCGT...ATATGCCAGCAGGCCGGCCACGCT
# outroot=data/processed/extract_seqs/barcode1

# ## Run cutadapt
# extract_seqs ${FQ} ${adapter} ${output}

# ## Barcode 2 ###################################################################
# adapter=ATATGCCAGCAGGCCGGCCACGCT...CGGTGGCCCGGGCGGCCGCACGATG
# outroot=data/processed/extract_seqs/barcode2

# ## Run cutadapt
# extract_seqs ${FQ} ${adapter} ${output}

# ## lacI ########################################################################
# adapter=TATTTTTTCCTCCTGGATTATCACT...TTCATATTCACCACCCTGAATTGAC
# outroot=data/processed/extract_seqs/lacI

# ## Run cutadapt
# extract_seqs ${FQ} ${adapter} ${output}

# ## Insulator ###################################################################
# adapter=CATCGTATAACGTTACTGGTTTCAT...GCTACATGAGTAGCAGTACGAAAAT
# outroot=data/processed/extract_seqs/insulator

# ## Run cutadapt
# extract_seqs ${FQ} ${adapter} ${output}

# ## tetA ########################################################################
# adapter=ACTAGCCATCAAGGAGAGCTGCTAC...GTGCCTAACGGCGTAAGGAGGTATT
# outroot=data/processed/extract_seqs/tetA

# ## Run cutadapt
# extract_seqs ${FQ} ${adapter} ${output}

# ## YEP #########################################################################
# adapter=CTAAACAAGAAACGAGTGCCTAACG...AATAAAAGCGGGAGACCAGAAACAA
# outroot=data/processed/extract_seqs/YEP

# ## Run cutadapt
# extract_seqs ${FQ} ${adapter} ${output}

# ## KAN #########################################################################
# adapter=TGGACGAGCTGTATAAATAAAAGCG...ACCCCTTAATAAGATGATCTTCTTG
# outroot=data/processed/extract_seqs/KAN

# ## Run cutadapt
# extract_seqs ${FQ} ${adapter} ${output}

# ## Ori #########################################################################
# adapter=TCCACTGAGCGTCAGACCCCTTAAT...ATAGTAAGCCAGTATACACTCCGCT
# outroot=data/processed/extract_seqs/Ori

# ## Run cutadapt
# extract_seqs ${FQ} ${adapter} ${output}

## Empty1 ######################################################################
adapter=ATCGGTGAGCCCGGGCTGTCGGCGT
outroot=data/processed/extract_seqs/empty1

## Run cutadapt
cutadapt \
    --rc --discard-untrimmed -j 14 \
    -a  ${adapter} \
    --info-file ${outroot}.info \
    -o ${outroot}.fastq ${FQ} \
    > ${outroot}.log

## Empty2 ######################################################################
adapter=CGGTGGCCCGGGCGGCCGCACGATG...TCACTGCCCGCTTTCCAGTCGGGAA 
outroot=data/processed/extract_seqs/empty2

## Run cutadapt
extract_seqs ${FQ} ${adapter} ${output}

## Empty3 ######################################################################
adapter=ACACCCTCATCAGTGCCAACATAGT
outroot=data/processed/extract_seqs/empty3

## Run cutadapt
cutadapt \
    --rc --discard-untrimmed -j 14 \
    --front ${adapter} \
    --info-file ${outroot}.info \
    -o ${outroot}.fastq ${FQ} \
    > ${outroot}.log

# ## leading_seq #################################################################
# adapter=ATATTTGCTCATGAGCCCGAAGTGG
# outroot=data/processed/extract_seqs/leading
# 
# ## Run cutadapt
# cutadapt \
#     --rc --discard-untrimmed -j 14 \
#     -a  ${adapter} \
#     --info-file ${outroot}.info \
#     -o ${outroot}.fastq ${FQ} \
#     > ${outroot}.log
# 
# ## trailing_seq ################################################################
# adapter=TGAGCGAGGAAGCACCTCAGATAAA
# outroot=data/processed/extract_seqs/trailing
# 
# ## Run cutadapt
# cutadapt \
#     --rc --discard-untrimmed -j 14 \
#     --front ${adapter} \
#     --info-file ${outroot}.info \
#     -o ${outroot}.fastq ${FQ} \
#     > ${outroot}.log

## NOTES ABOUT CUTADAPT OUTPUT
##### Info file - https://cutadapt.readthedocs.io/en/stable/guide.html#cutadapt-s-outputq
# The fields in a row that describes a match are:

# Read name
# Number of errors
# 0-based start coordinate of the adapter match
# 0-based end coordinate of the adapter match
# Sequence of the read to the left of the adapter match (can be empty)
# Sequence of the read that was matched to the adapter
# Sequence of the read to the right of the adapter match (can be empty)
# Name of the found adapter.
# Quality values corresponding to sequence left of the adapter match (can be empty)
# Quality values corresponding to sequence matched to the adapter (can be empty)
# Quality values corresponding to sequence to the right of the adapter match (can be empty)