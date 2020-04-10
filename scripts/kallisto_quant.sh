#!/bin/bash

index=genomes/Mmul_8.0.1.idx

# Assign a library ID
## samples.txt is a text file with each library ID listed (must be identical to prefixes used for fastq files)
lib=$(tail -n+2 data/rna_metadata.txt | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p output

# Map to transcriptome

## All reads are in the fastq/ directory with the naming convention ${lib}_S<sample_number>_L00<lane: 1 or 2>_R<read: 1 or 2>_001.fastq.gz

kallisto quant -i $index -t $SLURM_CPUS_ON_NODE -o output/$lib \
	fastq/$lib*L001*R1*.fastq.gz fastq/$lib_*L001*R2*.fastq.gz \
	fastq/$lib_*L002*R1*.fastq.gz fastq/$lib_*L002*R2*.fastq.gz
