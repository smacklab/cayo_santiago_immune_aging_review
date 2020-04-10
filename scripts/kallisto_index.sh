#!/bin/bash

mkdir -p genomes

# Fetch fasta from Ensembl
wget -O genomes/Mmul_8.0.1.fa.gz \
	ftp://ftp.ensembl.org/pub/release-96/fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_8.0.1.cdna.all.fa.gz

# Unzip fasta
gunzip genomes/Mmul_8.0.1.fa.gz

# Index fasta
kallisto index -i genomes/Mmul_8.0.1.idx genomes/Mmul_8.0.1.fa
