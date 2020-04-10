#!/usr/bin/env Rscript

source('scripts/_include_functions.R')

# Rename previously modeled rhesus effects
mmul.results = results

# Infidium HumanMethylation450 BeadChip metadata obtained from Illumina
##  (https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html)
## In particular, the following manifest file was used:
##   HumanMethylation450_15017482_v1-2.csv

# Data from Hannum et al. 2013 downloaded from the NCBI Gene Expression Omnibus (GSE40279)
## In particular, average beta values from the methylation array were used:
##   GSE40279_average_beta.txt.gz	

## Because of large file sizes and memory/time constraints, the files above can be
##   processed to remove excess columns and sites (i.e., sites that have no orthologs in
##   the comparative macaque data

# Import Infidium HumanMethylation450 BeadChip coordinates
hm = read.csv('data/HumanMethylation450_coords.csv',stringsAsFactors=FALSE)
rownames(hm) = hm[[1]]

# Import CpG IDs from the Hannum et al. 2013 dataset (GSE40279)
hannum.sites = scan(file='data/hannum_sites.txt',what='')

# Extract Infidium HumanMethylation450 coordinates for Hannum sites
hm.out = hm[hannum.sites,]

# Make site ID column
hm.out$bed = with(hm.out,paste0('chr',Chromosome_36,':',Coordinate_36,'-', Coordinate_36))
hm.out = subset(hm.out,Chromosome_36 != 'MULTI')

# Import orthologous sites in the rhesus dataset but in human coordinates (hg18)
rh = scan(what='',file='bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8_hg19_hg18.bed')

# Import the same sites but in rhesus coordinates (Mmul_8.0.1)
rh.rm8 = scan(what='',file='bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8.bed')

# Identify intersecting sites
sites.intersect = intersect(rh,hm.out$bed)

# Make data frame with Infidium sites and corresponding coordinates
df.hm = data.frame(id=hm.out$bed,ilmn=hm.out$IlmnID,stringsAsFactors=FALSE)

# Make data frame with rhesus sites and rhesus coordinates
df.rh = data.frame(id=rh,mmul=rh.rm8,stringsAsFactors=FALSE)

# Split chromosomes and position
df.hm = within(df.hm,{
	chr = gsub('chr(.+?):([0-9]+?)-([0-9]+)','\\1',id)
	pos = as.integer(gsub('chr(.+?):([0-9]+?)-([0-9]+)','\\2',id))
})

df.rh = within(df.rh,{
	chr = gsub('chr(.+?):([0-9]+?)-([0-9]+)','\\1',id)
	pos = as.integer(gsub('chr(.+?):([0-9]+?)-([0-9]+)','\\2',id))
})

df.rh = subset(df.rh,df.rh$chr %in% df.hm$chr)

# Load GenomicRanges to match corresponding sites and nearby sites

library(GenomicRanges)

i.rh = with(df.rh,GRanges(chr,IRanges(pos,width=1,names=id),'*'))
i.hm = with(df.hm,GRanges(chr,IRanges(pos,width=1,names=id),'*'))

olaps = nearest(i.rh, i.hm)

hits = data.frame(df.rh,
	chr.nearest = df.hm[olaps,]$chr,
	pos.nearest = df.hm[olaps,]$pos,
	ilmn.nearest = df.hm[olaps,]$ilmn,
	stringsAsFactors = FALSE)

# Identify nearest site
hits$id.nearest = with(hits,paste0('chr',chr.nearest,':',pos.nearest,'-',pos.nearest))

# Calculate distance to nearest site
hits$dist.nearest = with(hits,abs(pos-pos.nearest))

# Good hits are hits that are considered orthologous
## dist.nearest can be increased to include more sites (sites within 1kb tend to be highly comethylated)
good.hits = subset(hits,dist.nearest <= 100)

mmul.results$mmul8 = with(mmul.results,paste0('chr',chr,':',as.integer(pos),'-',as.integer(pos)))

good.hits.match = subset(good.hits,mmul %in% mmul.results$mmul8)

# Reduce macaque model results to essential columns and rename to reduce ambiguity
mmul.results = within(mmul.results,{
	beta.mmul = beta.age
	se.mmul = se.age
	pval.mmul = pval.age
	sbeta.mmul = sbeta.age
	fdr.mmul = fdr.age
})
mmul.results = subset(mmul.results,select=c('mmul8','beta.mmul','se.mmul','pval.mmul','sbeta.mmul','fdr.mmul'))

# Import average betas from human dataset (Hannum et al. 2013)
## hannum_good_hits.txt is a more manageable subset of the average betas file that contains matching sites to rhesus (based on a 100-bp neighboring site criterion)
hannum = read.delim('data/hannum_good_hits_100.txt',stringsAsFactors=FALSE)
hannum = hannum[complete.cases(hannum),]

# Model methylation is an analagous manner (see scripts/pqlseq) using a linear model
p.meth = as.matrix(hannum[2:ncol(hannum)])
rownames(p.meth) = hannum$ID_REF

keep = apply(p.meth,1,function(x) median(x) > 0.1 & median(x) < 0.9)

# Read in Hannum et al. samples metadata
meta = read.table('data/hannum_metadata.txt',sep='\t')
names(meta) = c('accession','id','age','source','plate','sex','race','tissue')

if (!all(colnames(p.meth) == meta$V2)) stop('Metadata do not match!')

lib.counts = matrix(1,nrow=nrow(p.meth),ncol=ncol(p.meth),dimnames=dimnames(p.meth))

library(parallel)
library(doParallel)

hsap.results = do.call(rbind,mclapply(1:nrow(p.meth),function(x) {
	d = data.frame(meta,int=p.meth[x,])
	model = lm(int~age + plate + sex + race,data=d)
	out = do.call(data.frame,as.list(as.numeric(summary(model)$coefficients[2:5,c(1,2,4)])))
	names(out) = c('beta.age','beta.plate','beta.sex','beta.race','se.age','se.plate','se.sex','se.race','pval.age','pval.plate','pval.sex','pval.race')
	rownames(out) = rownames(p.meth)[x]
	out
},mc.cores=8))

hsap.results$fdr.age = p.adjust(hsap.results$p.age,'fdr')

# Rename columns
hsap.results = within(hsap.results,{
	beta.hsap = beta.age
	se.hsap = se.age
	pval.hsap = p.age
	sbeta.hsap = beta.age / sqrt(se.age)
	fdr.hsap = p.adjust(pval.hsap,'fdr')
})
hsap.results = subset(hsap.results,select=c('beta.hsap','se.hsap','pval.hsap','sbeta.hsap','fdr.hsap'))

hsap.beta = hsap.results$beta.hsap
hsap.sbeta = hsap.results$sbeta.hsap
hsap.fdr = hsap.results$fdr.hsap

names(hsap.beta) = names(hsap.sbeta) = names(hsap.fdr) = rownames(hsap.results)

mmul.beta = results$beta_age
mmul.sbeta = with(results,beta_age / se_age)
mmul.fdr = results$fdr_age

names(mmul.beta) = names(mmul.sbeta) = names(mmul.fdr) = mmul.results$mmul8


good.hits.match$hsap.beta = hsap.beta[good.hits.match$ilmn.nearest]
good.hits.match$hsap.fdr = hsap.fdr[good.hits.match$ilmn.nearest]
good.hits.match$hsap.sbeta = hsap.sbeta[good.hits.match$ilmn.nearest]
good.hits.match$mmul.beta = mmul.beta[good.hits.match$mmul]
good.hits.match$mmul.fdr = mmul.fdr[good.hits.match$mmul]
good.hits.match$mmul.sbeta = mmul.sbeta[good.hits.match$mmul]

good.hits.match = good.hits.match[complete.cases(good.hits.match),]

bootstrap.concordance = function(cutoff) {
	concordance = with(subset(good.hits.match,mmul.fdr <= cutoff & hsap.fdr <= cutoff),(mmul.sbeta * hsap.sbeta) > 0)
	out = replicate(1000,mean(sample(concordance,replace=TRUE),na.rm=TRUE))
	data.frame(
		fdr = cutoff,
		concordant = sum(concordance,na.rm=TRUE),
		correlation = with(subset(good.hits.match,mmul.fdr <= cutoff & hsap.fdr <= cutoff),cor(mmul.beta,hsap.beta)),
		concordance = mean(concordance,na.rm=TRUE),
		stdev = sd(concordance,na.rm=TRUE),
		n = length(concordance),
		error.min = quantile(out,0.025,na.rm=TRUE),
		error.max = quantile(out,0.975,na.rm=TRUE),
		stringsAsFactors = FALSE)
}

hsap.mmul.concordance = do.call(rbind,lapply(seq(0,1,0.01),bootstrap.concordance))
