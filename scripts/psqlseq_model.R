#!/usr/bin/env Rscript

metadata = read.delim('data/dnam_metadata.txt',stringsAsFactors=FALSE)

# Convert dates
metadata = within(metadata,{
	dob = as.Date(dob)
	trap_date = as.Date(sample_date)
	age = as.numeric(sample_date - dob) / 365.25
	sex = factor(sex)
	batch = factor(batch)
})

p.meth = meth / total_counts

keep = apply(p.meth,1,function(x) median(x) > 0.1 & median(x) < 0.9)

meth = meth[keep,]
total_counts = total_counts[keep,]
p.meth = p.meth[keep,]

library(PQLseq)

dimnames(p.meth) = dimnames(meth) = dimnames(total_counts)

# Calculate kinship matrix from pedigree
bigped = read.delim('data/cayo_pedigree.txt',stringsAsFactors=FALSE)
names(bigped) = c('animal_id','sire','dam')

library(synbreed)
ped = create.pedigree(bigped$animal_id,bigped$sire,bigped$dam)
gp = create.gpData(pedigree=ped)
rmatrix = kin(gp,ret='add')
rmatrix = rmatrix[metadata$animal_id,metadata$animal_id]
dimnames(rmatrix) = list(metadata$library_id,metadata$library_id)
diag(rmatrix) = 1

# Set the resulting relatedness matrix to variable k
k = rmatrix

library(parallel)
library(doParallel)

results = pqlseq(
	RawCountDataSet = meth,
	Phenotypes = metadata$age,
	Covariates = model.matrix(~batch+sex,data=metadata),
	RelatednessMatrix = k,
	LibSize = total_counts,
	fit.model = 'BMM',
	filtering = TRUE,
	numCore=40
)

names(results) = gsub('_predictor$','_age',names(meth.model))

results = within(results,{
	beta.age = beta_age
	se.age = se_age
	pval.age = pval_age
	sbeta.age = beta / sqrt(se_age)
})

results$chr = gsub('(chr[0-9XY]+):([0-9]+)-([0-9]+)','\\1',rownames(results))
results$pos = gsub('(chr[0-9XY]+):([0-9]+)-([0-9]+)','\\2',rownames(results))

# Keep only converged and autosomal sites
results = subset(results,converged & !chr %in% c('chrX','chrY'))

dir.create('bed',showWarnings=FALSE,recursive=TRUE)

write(rownames(results),file='bed/rhesus_rrbs.bed',sep='\t')
