#!/usr/bin/bash
REF=$1
FQ=$2
OUTPREFIX=$3

minimap2 -a \
        -k 19 \
        -O 5,56 \
        -E 4,1 \
        -B 5 \
        -z 400,50 \
        -r 2k \
        --eqx \
        --secondary=no \
        -t 12 \
        ${REF} \
        ${FQ} \
    | samtools sort -@4 -m 10G -o ${OUTPREFIX}.bam 

# minimap2 -a \
#         -x splice -t 8 \
#         ${REF} \
#         ${FQ} |
#         samtools sort -@4 -m 10G -o ${OUTPREFIX}.bam 

samtools index ${OUTPREFIX}.bam