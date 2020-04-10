#!/usr/bin/env Rscript

# Fetch related GO terms
mmul.go = getBM(
	attributes=c('ensembl_gene_id', 'go_id', 'name_1006'),
	filters = 'ensembl_gene_id',
	values = rownames(results),
	mart = mmul)

# Assign related GO terms for each gene
gene2go = lapply(unique(mmul.go$ensembl_gene_id),function(x){
	out = sort(mmul.go[mmul.go$ensembl_gene_id==x,'go_id'])
	out[grepl('^GO:',out)]
})

names(gene2go)=unique(mmul.go$ensembl_gene_id)

# Save standardized betas
b = results$sbeta.age

names(b) = rownames(results)

library(topGO)

age.up.beta = new('topGOdata',
	description = 'Simple session',
	ontology = 'BP',
	allGenes = -b,
	geneSelectionFun = function(x) x > 0,
	nodeSize = 10,
	annotationFun = annFUN.gene2GO,
	gene2GO = gene2go)

age.up.go.ks = runTest(age.up.beta,algorithm='weight01',statistic='KS')

age.down.beta = new('topGOdata',
	description = 'Simple session',
	ontology = 'BP',
	allGenes = b,
	geneSelectionFun = function(x) x > 0,
	nodeSize = 10,
	annotationFun = annFUN.gene2GO,
	gene2GO = gene2go)

age.down.go.ks = runTest(age.down.beta,algorithm='weight01',statistic='KS')

# Create a join vector with names as the values and GO IDs as the indexes
go.join = mmul.go$name_1006
names(go.join) = mmul.go$go_id

library(reshape2)

age.b.up = melt(age.up.go.ks@score[p.adjust(age.up.go.ks@score,'fdr') < 0.1])
age.b.up$name = go.join[rownames(age.b.up)]
age.b.up$fdr = p.adjust(age.up.go.ks@score,'fdr')[p.adjust(age.up.go.ks@score,'fdr') < 0.1]

age.b.down = melt(age.down.go.ks@score[p.adjust(age.down.go.ks@score,'fdr') < 0.1])
age.b.down$name = go.join[rownames(age.b.down)]
age.b.down$fdr = p.adjust(age.down.go.ks@score,'fdr')[p.adjust(age.down.go.ks@score,'fdr') < 0.1]

## Some GO terms do not get matched to names.
##   The AMIGO website is more updated, so fill in these missing names using a web search
##   (parsing the HTML results as XML)

# Fetch missing GO names from the AMIGO site
to.do.b.up = which(is.na(age.b.up$name))
to.do.b.down = which(is.na(age.b.down$name))

library(tidyr)
library(XML)
for (i in to.do.b.up) {
	cat(i,'\n')
	search.url = paste0('http://amigo.geneontology.org/amigo/term/',rownames(age.b.up)[i])
	age.b.up$name[i] = search.url %>% htmlParse %>% xmlChildren %>% `[[`(3) %>% xmlChildren %>% `[[`(4) %>% xmlChildren %>% `[[`(24) %>% xmlChildren %>% `[[`(6) %>% xmlChildren %>% `[[`(2) %>% xmlChildren %>% `[[`('text') %>% xmlValue
}
for (i in to.do.b.down) {
	cat(i,'\n')
	search.url = paste0('http://amigo.geneontology.org/amigo/term/',rownames(age.b.down)[i])
	age.b.down$name[i] = search.url %>% htmlParse %>% xmlChildren %>% `[[`(3) %>% xmlChildren %>% `[[`(4) %>% xmlChildren %>% `[[`(24) %>% xmlChildren %>% `[[`(6) %>% xmlChildren %>% `[[`(2) %>% xmlChildren %>% `[[`('text') %>% xmlValue
}

# GO terms associated with increased expression
age.b.up = age.b.up[order(age.b.up$value),]

# GO terms associated with decreased expression
age.b.down = age.b.down[order(age.b.down$value),]
