# Rhesus macaque and human aging

This repository contains code written for a currently unpublished comparison of human and rhesus macaque aging effects in the peripheral blood immune system.

> Chiou KL, Montague MJ, Goldman EA, Watowich MM, Sams SN, Song J, Horvath JE, Sterner, KN, Ruiz-Lambides AV, Martínez MI, Higham JP, Brent LJN, Platt ML, Snyder-Mackler N. Rhesus macaques as a tractable physiological model of human ageing.

Citation details will be updated with the status of the paper.

Note that we ran most steps on the University of Washington high-performance computing cluster (Mox). We have aimed to generalize the code here by removing system-specific references to installed software and modules. Instead, we document required software and version numbers below (excluding standard Unix programs). For HPC systems, the required scripts and binaries must be in the PATH. The easiest way to do this is to use an existing module or to install your own. In these cases, the modules should be loaded prior to running the appropriate code below.

# Pipeline

## Gene expression analysis

1. First download and index the reference genome.

	Required software: [Kallisto](https://pachterlab.github.io/kallisto) (v0.43.1)

```
scripts/kallisto_index.sh
```

2. Next map reads and quantify expression.

	Required software: [Kallisto](https://pachterlab.github.io/kallisto) (v0.43.1)

```
sbatch --array=1-$(wc -l data/rna_samples.txt | cut -d ' ' -f 1) scripts/kallisto_quant.sh
```

3. Import expression matrix and model age effects.

	Required software: [R](https://cran.r-project.org) (v3.6.1)

```
scripts/kallisto_import.R
scripts/emma_model.R
```

4. Finally, perform enrichment analysis and comparison to human age effects.

	Required software: [R](https://cran.r-project.org) (v3.6.1)

```
scripts/topgo_enrichment.R
scripts/rna_human_comparison.R
```

## DNA methylation analysis

1. First download and index the reference genome.

	Required software: [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark) (v0.20.0)

```
scripts/bismark_index.sh
```

2. Next, map reads and extract methylated CpGs.

	Required software: [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark) (v0.20.0), [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore) (v0.4.5), [Cutadapt](https://github.com/marcelm/cutadapt) (v1.15), [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2) (v2.3.3.1), [Samtools](http://www.htslib.org) (v1.9)

```
scripts/bismark_extract.sh
```

3. Import methylation data and model age effects

```
scripts/bismark_import.R
scripts/pqlseq_model.R
```

4. Convert macaque coordinates to human (Hg18) coordinates, as used in the Infinium HumanMethylation450 BeadChip (Illumina)

	Required software: [liftOver](https://genome-store.ucsc.edu) (2018-03-24), [R](https://cran.r-project.org) (v3.6.1)

```
scripts/liftover_orthologs.sh
```

5. Finally, perform comparison to human age effects

	Required software: [R](https://cran.r-project.org) (v3.6.1)

```
scripts/dnam_human_comparison.R
```
