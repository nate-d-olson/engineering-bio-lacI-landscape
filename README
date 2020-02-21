# RAW READ QC
Fastq file initial QC using `fastqc -t 6 XTACK.ccs.fastq`, version FastQC v0.11.8
- Notes: ~1M total reads, with read lengths ranging from 101 bp to 21kb, 80K reads have read lengths around the expected 5.5kb


# READ ALIGNMENT
CCS reads aligned to full plasmid sequence using minimap2, see `align_reads.sh`.
Used minimap2 parameters from GIAB pacbio ccs 10Kb README

__Log from read alignment__ 
[M::mm_idx_gen::0.008*0.32] collected minimizers
[M::mm_idx_gen::0.008*0.48] sorted minimizers
[M::main::0.009*0.48] loaded/built the index for 1 target sequence(s)
[M::mm_mapopt_update::0.009*0.49] mid_occ = 3
[M::mm_idx_stat] kmer size: 19; skip: 10; is_hpc: 0; #seq: 1
[M::mm_idx_stat::0.009*0.49] distinct minimizers: 1014 (99.90% are singletons); average occurrences: 1.001; average spacing: 5.577
[M::worker_pipeline::42.385*3.00] mapped 97347 sequences
[M::worker_pipeline::74.974*3.03] mapped 94453 sequences
[M::worker_pipeline::107.589*3.05] mapped 93946 sequences
[M::worker_pipeline::140.214*3.05] mapped 94077 sequences
[M::worker_pipeline::173.073*3.06] mapped 94041 sequences
[M::worker_pipeline::205.467*3.06] mapped 94131 sequences
[M::worker_pipeline::237.985*3.06] mapped 94522 sequences
[M::worker_pipeline::270.428*3.07] mapped 94423 sequences
[M::worker_pipeline::302.777*3.07] mapped 95123 sequences
[M::worker_pipeline::335.196*3.07] mapped 97015 sequences
[M::worker_pipeline::346.141*3.06] mapped 59121 sequences
[M::main] Version: 2.17-r941
[M::main] CMD: minimap2 -a -k 19 -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k --eqx --secondary=no XTACK_full_plasmid.fa XTACK.ccs.fastq
[M::main] Real time: 346.143 sec; CPU: 1058.936 sec; Peak RSS: 2.429 GB
[bam_sort_core] merging from 10 files and 1 in-memory blocks...

# Alignment QC
calculating bam stats using `samtools stats -r XTACK_full_plasmid.fa XTACK.ccs.bam > XTACK.ccs.bam.stats`

# Annotated plasmid
Aligned annotated fasta to full plasmid sequence
`minimap2 -a XTACK_full_plasmid.fa XTACK_annotated_plasmid.fa > XTACK_annotated.sam` 
Generated bed file from alignment
`bedtools bamtobed -i XTACK_annotated.sam > XTACK_annotated.bed`

# Reordering plasmid
Reordered plasmid sequence to start with tet as most a lot of the reads were clipped at the left end of the alignment. 
`XTACK_full_plasmid_reordered.fa`
`XTACK_annotated_reordered.bed`
Hoping to get better coverage of the barcode and LacI sequence. 

Reordering did not really help. Need to look into better method for read sequence alignment


## Extracting target sequences
__Approach__ 
Using cutadapt (https://cutadapt.readthedocs.io/en/stable/guide.html) to pull out targets

Target   | 5prime                    | 3prime                    |
barcode1 | ATCGGTGAGCCCGGGCTGTCGGCGT | ATATGCCAGCAGGCCGGCCACGCT  |
barcode2 | ATATGCCAGCAGGCCGGCCACGCT  | CGGTGGCCCGGGCGGCCGCACGATG |
lacI     | TATTTTTTCCTCCTGGATTATCACT | TTCATATTCACCACCCTGAATTGAC |
Insulator| CATCGTATAACGTTACTGGTTTCAT | GCTACATGAGTAGCAGTACGAAAAT |
tetA     | ACTAGCCATCAAGGAGAGCTGCTAC | GTGCCTAACGGCGTAAGGAGGTATT |
YEP      | CTAAACAAGAAACGAGTGCCTAACG | AATAAAAGCGGGAGACCAGAAACAA |
KAN      | TGGACGAGCTGTATAAATAAAAGCG | ACCCCTTAATAAGATGATCTTCTTG |
Ori      | TCCACTGAGCGTCAGACCCCTTAAT | ATAGTAAGCCAGTATACACTCCGCT |
empty1   | ATCGGTGAGCCCGGGCTGTCGGCGT ||
empty2   | CGGTGGCCCGGGCGGCCGCACGATG | TCACTGCCCGCTTTCCAGTCGGGAA |
empty3   | ACACCCTCATCAGTGCCAACATAGT ||
leading  | ATATTTGCTCATGAGCCCGAAGTGG ||
training | TGAGCGAGGAAGCACCTCAGATAAA ||


# Notes
- Might want to look into read mapper for circular genomes - or find out how to best deal with supplementary alignments
- filtering reads by length http://www.metagenomics.wiki/tools/short-read/filter-remove-too-short-reads-length