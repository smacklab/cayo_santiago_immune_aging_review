#!/bin/bash

mkdir -p genomes

# Fetch fasta from Ensembl
wget -O genomes/bismark/Mmul_8.0.1.fa.gz \
	ftp://ftp.ensembl.org/pub/release-96/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_8.0.1.dna.toplevel.fa.gz

# Unzip fasta
gunzip genomes/bismark/Mmul_8.0.1.fa.gz

# Index fasta
bismark_genome_preparation genomes/bismark
