#!/usr/bin/bash
PLASMIDSEQ=$1
ANNOPLASMIDSEQ=$2
OUTPREFIX=$3

minimap2 -a ${PLASMIDSEQ} ${ANNOPLASMIDSEQ} \
    | bedtools bamtobed -i - \
    | bedtools sort -i - \
    > $OUTPREFIX.bed