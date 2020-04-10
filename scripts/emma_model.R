#!/usr/bin/env Rscript

# Load in trapping/animal metadata
metadata = read.delim('data/rna_metadata.txt',stringsAsFactors=FALSE)

# Convert dates
metadata = within(metadata,{
	dob = as.Date(dob)
	trap_date = as.Date(sample_date)
	age = as.numeric(sample_date - dob) / 365.25
	sex = factor(sex)
	batch = factor(batch)
})

library(limma)
library(edgeR)

# Vector of genes with a TPM >= 2
keep = names(which(rowMeans(txi.gene$abundance)>=2))

# Normalize count matrix
v = voom(txi.gene$counts)

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

# Create z matrix as identity matrix

## Set dimensions and dimension names
z = matrix(0,nrow=nrow(metadata),ncol=nrow(k))
rownames(z) = metadata$library_id
colnames(z) = rownames(k)

## Set diagonal equal to 1
i = with(metadata,cbind(library_id,library_id))
z[i] = 1

# Set expression matrix passing filter
e = v$E[keep,]

library(parallel)
library(doParallel)

m = metadata

# Set covariates
design = model.matrix(~sex + age + batch,data=m,)

clus = makeCluster(4)
registerDoParallel(cores=4)
clusterExport(clus,varlist=c('e','k','z','design'),envir=environment())

results = t(parApply(clus,e,1,function(y) {
	require(EMMREML)

	# Run model
	emma=emmreml(y = y,X = design,Z = z,K = k,varbetahat = T,varuhat = T,PEVuhat = T,test = T)

	# Extract p value
	p = emma$pvalbeta

	# Extract variance of model beta
	varb = emma$varbetahat

	# Extract model beta
	b = emma$betahat

	c(b,varb,p[,"none"])
}))

stopCluster(clus)

# Name output
colnames(results)[(ncol(design) * 0 + 1):(ncol(design) * 1)] = paste('beta',colnames(design),sep='.')
colnames(results)[(ncol(design) * 1 + 1):(ncol(design) * 2)] = paste('bvar',colnames(design),sep='.')
colnames(results)[(ncol(design) * 2 + 1):(ncol(design) * 3)] = paste('pval',colnames(design),sep='.')

results = as.data.frame(results)

# Set standardized beta
results$sbeta.age = results$beta.age / sqrt(results$bvar.age)

# Set FDR (q value)
results$fdr.age = p.adjust(results$pval.age,'fdr')
