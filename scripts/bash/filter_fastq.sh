#!/usr/bin/bash

## Filter fastq by read length
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 3000 && length(seq) <= 6500) {print header, seq, qheader, qseq}}' $1
