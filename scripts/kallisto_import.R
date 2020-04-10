#!/usr/bin/env Rscript

library(tximport)
library(rhdf5)

samples = read.delim('data/rna_metadata.txt',stringsAsFactors=FALSE)[[1]]

files = file.path('output',samples,'abundance.h5')
names(files) = samples

## Import hd5 kallisto-mapped files
txi.kallisto = tximport(files, type = 'kallisto', txOut = TRUE)

library(biomaRt)
mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl')

# Fetch transcript-to-gene join table
tx2gene = getBM(
	attributes = c('ensembl_transcript_id_version','ensembl_gene_id'),
	filters = 'ensembl_transcript_id_version',
	values = rownames(txi.kallisto$abundance),
	mart = mmul)

# Summarize mapped data to gene level
txi.gene = summarizeToGene(txi.kallisto, tx2gene)
