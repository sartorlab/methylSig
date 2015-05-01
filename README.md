# methylSig

A whole genome DNA methylation analysis pipeline.

## Description

MethylSig is a method for testing for differentially methylated cytosines (DMCs) or regions (DMRs) in whole-genome bisulfite sequencing (bis-seq) or reduced representation bisulfite sequencing (RRBS) experiments. MethylSig uses a beta binomial model to test for significant differences between groups of samples. Several options exist for either site-specific or sliding window tests, combining strands, and for variance estimation. It allows annotating the resulting regions to multiple genome features, and visualizing the results for chosen genomic regions along with supporting genomic information.

## Installation

MethylSig is not currently on CRAN or Bioconductor. Installation is easiest with [devtools](http://cran.r-project.org/web/packages/devtools/index.html).

```R
library(devtools)
install_github(‘rcavalcante/methylSig')
```
