
# methylSig

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/sartorlab/methylSig.svg?branch=master)](https://travis-ci.org/sartorlab/methylSig)
[![Coveralls test coverage](https://coveralls.io/repos/github/sartorlab/methylSig/badge.svg)](https://coveralls.io/r/sartorlab/methylSig?branch=master)
<!-- badges: end -->

MethylSig is a method for testing for differentially methylated cytosines (DMCs) or regions (DMRs) in whole-genome bisulfite sequencing (WGBS) or reduced representation bisulfite sequencing (RRBS) experiments.  MethylSig uses a beta binomial model to test for significant differences between groups of samples. Several options exist for either site-specific or sliding window tests, and variance estimation.

## Installation

You can install the released version of methylSig from [GitHub](https://github.com/sartorlab/methylSig) with:

``` r
devtools::install_github("sartorlab/methylSig@Bioc_3_10")
```
