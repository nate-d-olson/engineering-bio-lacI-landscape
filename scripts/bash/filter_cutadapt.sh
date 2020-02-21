#!/usr/bin/bash

for fq in *fastq; do
    ## Filter fastq by read length
    awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 1 && length(seq) <= 6500) {print header, seq, qheader, qseq}}' $fq \
    > ${fq%.*}_filtered.fq
done