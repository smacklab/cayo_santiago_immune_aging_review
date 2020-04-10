#!/usr/bin/env Rscript

source('scripts/_include_functions.R')

# Fetch macaque-human homologs
mmul.hsap = getBM(
	attributes=c('ensembl_gene_id', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name'),
	filters = 'ensembl_gene_id',
	values = rownames(results),
	mart = mmul)

# Copy model output to new data frame and combine it with human metadata
hsap.results = results
hsap.results$ensembl_gene_id = rownames(hsap.results)
hsap.results = merge(hsap.results,mmul.hsap,by='ensembl_gene_id')

# To avoid confusion in column names, specify that generically named columns are for macaques
hsap.results = within(hsap.results,{
	beta.mmul = beta.age
	p.mmul = pval.age
	fdr.mmul = fdr.age
	sbeta.mmul = sbeta.age
})

# Import data from Peters et al. 2015 (doi: 10.1038/ncomms9570)
## The file can be downloaded as Supplementary Data 1 from the Peters et al. 2015 paper
##   and moved to the data/ folder

library(gdata)
peters = read.xls('data/41467_2015_BFncomms9570_MOESM436_ESM.xlsx',sheet=1,skip=2,header=TRUE)

# Select and rename relevant columns
peters = peters[c(5,6,9,11,14,16,18,20)]
names(peters) = c('entrez','gene_id','disc.z','disc.p','repl.z','repl.p','meta.z','meta.p')

peters$disc.fdr = p.adjust(peters$disc.p,'fdr')
peters$repl.fdr = p.adjust(peters$repl.p,'fdr')
peters$meta.fdr = p.adjust(peters$meta.p,'fdr')

# Filter to only genes with macaque homologs
peters = subset(peters,gene_id %in% mmul.hsap$hsapiens_homolog_associated_gene_name & !is.na(peters$meta.p))

# Use the metaanalysis statistics
peters = within(peters,{
	z.hsap = meta.z
	p.hsap = meta.p
	fdr.hsap = meta.fdr
})

# Join associated human and macaque genes
hsap.mmul = merge(hsap.results,peters,by.x='hsapiens_homolog_associated_gene_name',by.y='gene_id')

# Keep table passing FDR < 0.2 threshold in both datasets
hsap.mmul.sig = subset(hsap.mmul,fdr.mmul <= 0.2 & fdr.hsap <= 0.2)


#  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #
#                                                                       #
#    Analysis of directional concordance between macaques and humans    #
#                                                                       #
#  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #

# Concordance refers to the consistency in directional age effects
hsap.mmul$concordance = with(hsap.mmul,sbeta.mmul * z.hsap > 0)

# Perform 1000 bootstrap replicates across a range of FDR thresholds and summarize data
hsap.mmul.concordance = do.call(rbind,lapply(seq(0.01,1,0.01),function(x) {
	out = replicate(1000,mean(sample(subset(hsap.mmul,fdr.mmul<=x & fdr.hsap<=x)$concordance,replace=TRUE)))
	pass.filter = subset(hsap.mmul,fdr.mmul<=x & fdr.hsap<=x)
	data.frame(
		fdr = x,
		concordant = with(pass.filter,sum(concordance)),
		concordance = with(pass.filter,mean(concordance)),
		stdev = with(pass.filter,sd(concordance)),
		n = nrow(pass.filter),
		error.min = quantile(out,0.025),
		error.max = quantile(out,0.975),
		stringsAsFactors = FALSE)
}))


#  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #
#                                                                       #
#          Macaque age predictions from human blood predictors          #
#                                                                       #
#  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  #

# Center and scale the expression data
e.scaled = t(apply(e,1,scale))
dimnames(e.scaled) = list(rownames(e),colnames(e))

# Read in predictors from Peters et al. (doi: 10.1038/ncomms9570)
## The file can be downloaded as Supplementary Data 5 from the Peters et al. 2015 paper
##   and moved to the data/ folder
peters.predict = read.xls('data/41467_2015_BFncomms9570_MOESM440_ESM.xlsx',sheet=2,skip=2,header=TRUE)
peters.predict = peters.predict[c(3,4)]
names(peters.predict) = c('gene_id','z')

# Merge with macaque data and keep only shared genes
peters.predict = merge(peters.predict,mmul.hsap,by.x='gene_id',by.y='hsapiens_homolog_associated_gene_name')
peters.predict = subset(peters.predict,ensembl_gene_id %in% rownames(e.scaled))

# Get rid of genes that have multiple macaque ENSEMBL IDs
e.scaled = e.scaled[rownames(e.scaled) %in% hsap.mmul$ensembl_gene_id,]
peters.predict = subset(peters.predict,ensembl_gene_id %in% rownames(e.scaled))
peters.predict = rbind(
	subset(peters.predict,!ensembl_gene_id %in% names(which(table(peters.predict$ensembl_gene_id)>1))),
	do.call(rbind,
		lapply(split(subset(peters.predict,ensembl_gene_id %in% names(which(table(peters.predict$ensembl_gene_id)>1))),subset(peters.predict,ensembl_gene_id %in% names(which(table(peters.predict$ensembl_gene_id)>1)))$ensembl_gene_id),function(x) x[1,]))
)
rownames(peters.predict) = peters.predict$ensembl_gene_id

# Equalize gene orders
e.scaled = e.scaled[sort(rownames(e.scaled)),]
peters.predict = peters.predict[rownames(e.scaled),]

# Calculate predicted ages (Z) as the product of the predictor and the normalized read count for each gene
age.z = unlist(lapply(colnames(e.scaled),function(x) sum(e.scaled[,x] * peters.predict$z)))
names(age.z) = colnames(e.scaled)

# Calculated standardized age predictions (Zs) to scale to actual rhesus ages
age.zs = mean(metadata$age) + (age.z - mean(age.z)) * sd(metadata$age) / sd(age.z)

rhesus.predictions = subset(metadata,select=c('animal_id','library_id','sex','age'))

# Redo predictions across a range of FDR thresholds
rhesus.predictions.thresholds = do.call(rbind,lapply(seq(0.01,1,0.01),function(x) {
	e.scaled.fdr = e.scaled[rownames(e.scaled) %in% subset(hsap.mmul,fdr.mmul <= x & fdr.hsap <= x)$ensembl_gene_id,]
	peters.predict.fdr = peters.predict[rownames(e.scaled.fdr),]
	age.z.fdr = unlist(lapply(colnames(e.scaled.fdr),function(x) sum(e.scaled.fdr[,x] * peters.predict.fdr$z)))
	names(age.z.fdr) = colnames(e.scaled.fdr)
	age.zs.fdr = mean(metadata$age) + (age.z.fdr - mean(age.z.fdr)) * sd(metadata$age) / sd(age.z.fdr)
	out = rhesus.predictions
	out$fdr = x
	out$zs = age.zs.fdr[out$library_id]
	out
}))
