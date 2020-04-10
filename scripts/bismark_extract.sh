#!/bin/bash

lib=$(tail -n+2 data/dnam_metadata.txt | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p trimmed

# Trim adaptors using Trim Galore! and RRBS specific parameters
trim_galore --paired --rrbs --non_directional --gzip \
	-o trimmed fastq/${lib}.R1.fastq.gz fastq/${lib}.R2.fastq.gz

mkdir -p bam

bismark --genome_folder genomes/bismark -o bam/ --non_directional --parallel 8 --score_min L,-0.6,-0.6 -1 trimmed/${lib}.R1_val_1.fq.gz -2 trimmed/${lib}.R2_val_1.fq.gz

mkdir -p output
mkdir -p output/Mmul_8.0.1

bismark_methylation_extractor -p --no_overlap --bedGraph --multicore 2 \
	--comprehensive --merge_non_CpG -o output/bismark \
	--genome_folder genomes/bismark bam/${lib}.R1_val_1_bismark_bt2_pe.bam
