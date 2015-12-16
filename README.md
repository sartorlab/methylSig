# methylSig

[![Build Status](https://travis-ci.org/sartorlab/methylSig.svg?branch=master)](https://travis-ci.org/sartorlab/methylSig) [![Coverage Status](https://coveralls.io/repos/sartorlab/methylSig/badge.svg?branch=master&service=github)](https://coveralls.io/github/sartorlab/methylSig?branch=master)

A whole genome DNA methylation analysis pipeline.

## Description

MethylSig is a method for testing for differentially methylated cytosines (DMCs) or regions (DMRs) in whole-genome bisulfite sequencing (bis-seq) or reduced representation bisulfite sequencing (RRBS) experiments. MethylSig uses a beta binomial model to test for significant differences between groups of samples. Several options exist for either site-specific or sliding window tests, combining strands, and for variance estimation. It allows annotating the resulting regions to multiple genome features, and visualizing the results for chosen genomic regions along with supporting genomic information.

## Installation

MethylSig is not currently on CRAN or Bioconductor. Installation is easiest with [devtools](http://cran.r-project.org/web/packages/devtools/index.html).

```R
library(devtools)
install_github('sartorlab/methylSig')
```

## Citation

Yongseok Park, Maria E. Figueroa, Laura S. Rozek, and Maureen A. Sartor, MethylSig: a whole genome DNA methylation analysis pipeline, *Bioinformatics* (2014) 30 (17): 2414-2422, doi:[10.1093/bioinformatics/btu339](http://bioinformatics.oxfordjournals.org/content/30/17/2414)
