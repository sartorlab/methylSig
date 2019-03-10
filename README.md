# methylSig

[![Build Status](https://travis-ci.org/sartorlab/methylSig.svg?branch=master)](https://travis-ci.org/sartorlab/methylSig) [![Coverage Status](https://coveralls.io/repos/github/sartorlab/methylSig/badge.svg?branch=master)](https://coveralls.io/github/sartorlab/methylSig?branch=master)

A whole genome DNA methylation analysis pipeline.

## Description

MethylSig is a method for testing for differentially methylated cytosines (DMCs) or regions (DMRs) in whole-genome bisulfite sequencing (bis-seq) or reduced representation bisulfite sequencing (RRBS) experiments. MethylSig uses a beta binomial model to test for significant differences between groups of samples. Several options exist for either site-specific or sliding window tests, combining strands, and for variance estimation. It allows annotating the resulting regions to multiple genome features, and visualizing the results for chosen genomic regions along with supporting genomic information.

## Installation

`methylSig` is available on GitHub at <http://www.github.com/sartorlab/methylSig>. Please use the appropriate branch for your version of R / Bioconductor:

## R 3.4 and Bioconductor 3.6

```{r, eval=FALSE}
devtools::install_github('sartorlab/methylSig@R3.4_Bioc3.6')
```

## R 3.5 and Bioconductor 3.7

Forthcoming.

## R 3.5 and Bioconductor 3.8

```{r, eval=FALSE}
devtools::install_github('sartorlab/methylSig@R3.5_Bioc3.8')
```

## Citation

Yongseok Park, Maria E. Figueroa, Laura S. Rozek, and Maureen A. Sartor, MethylSig: a whole genome DNA methylation analysis pipeline, *Bioinformatics* (2014) 30 (17): 2414-2422, doi:[10.1093/bioinformatics/btu339](http://bioinformatics.oxfordjournals.org/content/30/17/2414)
