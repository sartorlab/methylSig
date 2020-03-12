#' MethylSig: Differential Methylation Testing for WGBS and RRBS Data
#'
#' MethylSig is a package for testing for
#'    differentially methylated cytosines (DMCs) or regions (DMRs) in
#'    whole-genome bisulfite sequencing (WGBS) or reduced representation
#'    bisulfite sequencing (RRBS) experiments.  MethylSig uses a beta
#'    binomial model to test for significant differences between groups of
#'    samples. Several options exist for either site-specific or sliding
#'    window tests, and variance estimation.
#'
#' @section methylSig functions:
#' filter_loci_by_coverage()
#' filter_loci_by_snps()
#' tile_by_regions()
#' tile_by_windows()
#' filter_loci_by_group_coverage()
#' diff_binomial()
#' diff_methylsig()
#' diff_methylsig_dss()
#' annotate_diff()
#' visualize_diff()
#' region_enrichment_diff()
#'
#' @import bsseq
#' @import DelayedArray
#' @import DelayedMatrixStats
#' @import DSS
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import IRanges
#' @import parallel
#' @importFrom methods is
#' @importFrom stats as.formula
#' @importFrom stats formula
#' @importFrom stats pchisq
#' @importFrom stats pt
#' @importFrom stats p.adjust
#'
#' @docType package
#' @name methylSig
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
