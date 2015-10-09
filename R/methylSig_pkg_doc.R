#' methylSig: a whole genome DNA methylation analysis pipeline
#'
#' MethylSig is a method for testing differential methylated cytosines (DMCs) or regions (DMRs) in whole-genome bisulfite sequencing (bis-seq) or reduced representation bisulfite sequencing (RRBS) experiments. MethylSig uses a beta binomial model to test for significant differences between groups of samples. Several options exist for either site-specific or sliding window tests, combining strands, and for variance estimation. It allows annotating the resulting regions to multiple genome features, and visualizing the results for chosen genomic regions.
#'
#' @author Yongseok Park \email{yongpark@@pitt.edu}, Raymond Cavalcante \email{rcavalca@@umich.edu}, and Maureen A. Sartor
#' @references \url{http://www.github.com/sartorlab/methylSig}{MethylSig}
#'
#' @importFrom parallel mclapply
#' @importFrom boot corr
#' @import GenomicRanges
#' @import methods
#'
#' @docType package
#' @name methylSig-package
#' @aliases methylSig
NULL
