#!/usr/bin/bash

REF=$1
FQ=$2
OUTPREFIX=$3

graphmap align --min-read-len 1000 \
    -t 12 -C -r ${REF} -d ${FQ} -o ${OUTPREFIX}.sam

samtools sort -@4 -m 10G -o ${OUTPREFIX}.bam ${OUTPREFIX}.sam

samtools index ${OUTPREFIX}.bam

