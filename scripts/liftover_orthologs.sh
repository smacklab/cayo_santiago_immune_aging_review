#!/bin/bash

mkdir -p chain

# Fetch chain files from UCSC Genome Browser
wget -O chain/rheMac8ToHg19.over.chain.gz https://hgdownload-test.gi.ucsc.edu/goldenPath/rheMac8/liftOver/rheMac8ToHg19.over.chain.gz
wget -O chain/hg19ToRheMac8.over.chain.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToRheMac8.over.chain.gz
wget -O chain/hg19ToHg18.over.chain.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg18.over.chain.gz
wget -O chain/hg18ToHg19.over.chain.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz

for i in `ls chain/*.gz`; do gunzip $i; done

# Convert from Mmul8 to hg19
liftOver -positions bed/rhesus_rrbs.bed chain/rheMac8ToHg19.over.chain bed/rhesus_rrbs_hg19.mapped bed/rhesus_rrbs_hg19.unmapped
scripts/rmdup bed/rhesus_rrbs_hg19.mapped

# Back-convert to Mmul8
liftOver -positions rhesus_rrbs_hg19.bed chain/hg19ToRheMac8.over.chain bed/rhesus_rrbs_hg19_rm8.mapped bed/rhesus_rrbs_hg19_rm8.unmapped
scripts/rmdup bed/rhesus_rrbs_hg19_rm8.mapped

# Convert to hg19
liftOver -positions rhesus_rrbs_hg19_rm8.bed chain/rheMac8ToHg19.over.chain bed/rhesus_rrbs_hg19_rm8_hg19.mapped bed/rhesus_rrbs_hg19_rm8_hg19.unmapped
scripts/rmdup bed/rhesus_rrbs_hg19_rm8_hg19.mapped

# Convert to hg18
liftOver -positions rhesus_rrbs_hg19_rm8_hg19.bed chain/hg19ToHg18.over.chain bed/rhesus_rrbs_hg19_rm8_hg19_hg18.mapped bed/rhesus_rrbs_hg19_rm8_hg19_hg18.unmapped
scripts/rmdup bed/rhesus_rrbs_hg19_rm8_hg19_hg18.mapped

# Back-convert to hg19
liftOver -positions rhesus_rrbs_hg19_rm8_hg19_hg18.bed chain/hg18ToHg19.over.chain bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19.mapped bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19.unmapped
scripts/rmdup bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19.mapped

# Convert to hg18
liftOver -positions rhesus_rrbs_hg19_rm8_hg19_hg18_hg19.bed chain/hg19ToHg18.over.chain bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18.mapped bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18.unmapped
scripts/rmdup bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18.mapped

# Convert to hg19
liftOver -positions rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18.bed chain/hg18ToHg19.over.chain bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19.mapped bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19.unmapped
scripts/rmdup bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19.mapped

# Convert to Mmul8
liftOver -positions rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19.bed chain/hg19ToRheMac8.over.chain bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8.mapped bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8.unmapped
scripts/rmdup bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8.mapped

# At this point, only single-copy orthologs should remain that are common to all three assemblies.
# To finish, generate final versions in coordinates of hg19 and especially hg18

# Convert to hg19
liftOver -positions rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8.bed chain/rheMac8ToHg19.over.chain bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8_hg19.mapped bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8_hg19.unmapped
scripts/rmdup bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8_hg19.mapped

# Convert to hg18
liftOver -positions rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8_hg19.bed chain/hg19ToHg18.over.chain bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8_hg19_hg18.mapped bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8_hg19_hg18.unmapped
scripts/rmdup bed/rhesus_rrbs_hg19_rm8_hg19_hg18_hg19_hg18_hg19_rm8_hg19_hg18.mapped
