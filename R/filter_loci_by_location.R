#' Remove loci by overlap with a \code{GRanges} object
#'
#' A function to remove loci from a \code{BSseq} object based on intersection with loci in a \code{GRanges} object.
#'
#' @param bs a \code{BSseq} object.
#' @param gr a \code{GRanges} object.
#'
#' @return A \code{BSseq} object with loci intersecting \code{gr} removed.
#'
#' @examples
#' data(bsseq_stranded, package = 'methylSig')
#' regions = GenomicRanges::GRanges(
#'     seqnames = c('chr1','chr1','chr1','chr1'),
#'     ranges = IRanges::IRanges(
#'         start = c(5,25,45,70),
#'         end = c(15,40,55,80)
#'     )
#' )
#' filtered = filter_loci_by_location(bs = bsseq_stranded, gr = regions)
#'
#' @export
filter_loci_by_location = function(bs, gr) {

    # Missing checks
    if (missing(bs)) {
        stop('Must pass bs as a BSseq object.')
    }
    if (missing(gr)) {
        stop('Must pass gr as a GRanges object.')
    }

    # Type checks
    if (!is(bs, 'BSseq')) {
        stop('bs must be class BSseq.')
    }
    if (!is(gr, 'GRanges')) {
        stop('gr must be class GRanges.')
    }

    #####################################

    overlaps = GenomicRanges::findOverlaps(bs, gr)

    keep_idx = setdiff(seq_along(bs), unique(S4Vectors::queryHits(overlaps)))

    if(length(keep_idx) == 0) {
        stop('All loci in bs were removed by gr, leaving no loci for downstream analysis.')
    }

    bs = bs[keep_idx]

    return(bs)
}
