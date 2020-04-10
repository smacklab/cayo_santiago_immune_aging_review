#!/usr/bin/env Rscript

samples = read.delim('data/dnam_metadata.txt',stringsAsFactors=FALSE)[[1]]

library(parallel)
options(stringsAsFactors=FALSE)

chr = c(1:21,'X','Y')

cpg.bismark = mclapply(samples,function(x) {
	out = read.table(gzfile(file.path('output/bismark',paste0(x,'.R1_val_1_bismark_bt2_pe.bismark.cov.gz'))),header=FALSE,skip=1)
	out = subset(out,V1 %in% chr)
	out$chr = as.character(out$V1)
	out$id = with(out,paste0('chr',chr,':',V2,'-',V2))
	out$mcpg = with(out,V5)
	out$reads = with(out,V5 + V6)
	result = cbind(out$mcpg,out$reads)
	rownames(result) = out$id
	colnames(result) = c('mcpg','reads')
	result
},mc.cores=40)

all_sites = unique(unlist(lapply(total_reads,rownames)))

meth = do.call(cbind,lapply(cpg.bismark,function(x) {
	out = integer(length(poo))
	names(out) = all_sites
	out[rownames(x)] = x[,'mcpg']
	out
}))

total_counts = do.call(cbind,lapply(cpg.bismark,function(x) {
	out = integer(length(poo))
	names(out) = all_sites
	out[rownames(x)] = x[,'reads']
	out
}))

dimnames(meth) = dimnames(total_counts) = list(all_sites,samples)
