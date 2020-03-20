#' BSseq object read from stranded coverage files
#'
#' Data contains 11 methylation loci and 2 samples
#'
#' @format A BSseq object
#' @source data-raw/02-create_bsseq_rda.R
#' @examples
#' data(bsseq_stranded, package = 'methylSig')
'bsseq_stranded'

#' BSseq object read from destranded coverage files
#'
#' Data contains 6 methylation loci and 2 samples
#'
#' @format A BSseq object
#' @source data-raw/02-create_bsseq_rda.R
#' @examples
#' data(bsseq_destranded, package = 'methylSig')
'bsseq_destranded'

#' BSseq object with loci on multiple chromosomes
#'
#' Data contains 4 methylation loci for 2 samples on 2 chromosomes
#'
#' @format A BSseq object
#' @source data-raw/02-create_bsseq_rda.R
#' @examples
#' data(bsseq_multichrom, package = 'methylSig')
'bsseq_multichrom'

#' GRanges object with collapsed promoters on chr21 and chr22
#'
#' Data contains 1466 promoters for use in the vignette
#'
#' @format A GRanges object
#' @source data-raw/02-create_bsseq_rda.R
#' @examples
#' data(promoters_gr, package = 'methylSig')
'promoters_gr'
