#' Filter BSseq object by coverage
#'
#' Used after \code{bsseq::read.bismark} to mark loci in samples below \code{min_count} or above \code{max_count} to 0. These loci will then be removed prior to differential analysis by \code{filter_loci_by_group_coverage()} if there are not a sufficient number of samples with appropriate coverage.
#'
#' @param bs a \code{BSseq} object resulting from \code{bsseq::read.bismark} or constructed manually by the user.
#' @param min_count an \code{integer} giving the minimum coverage required at a locus.
#' @param max_count an \code{integer} giving the maximum coverage allowed at a locus.
#'
#' @return A \code{BSseq} object with samples/loci in the coverage and methylation matrix set to 0 where the coverage was less than \code{min_count} or greater than \code{max_count}. The number of samples and loci are conserved.
#'
#' @examples
#' bis_cov_file1 = system.file('extdata', 'bis_cov1.cov', package = 'methylSig')
#' bis_cov_file2 = system.file('extdata', 'bis_cov2.cov', package = 'methylSig')
#' test = bsseq::read.bismark(
#'     files = c(bis_cov_file1, bis_cov_file2),
#'     colData = data.frame(row.names = c('test1','test2')),
#'     rmZeroCov = FALSE,
#'     strandCollapse = FALSE
#' )
#' test = filter_loci_by_coverage(bs = test, min_count = 10, max_count = 500)
#' @export
filter_loci_by_coverage = function(bs, min_count = 5, max_count = 500) {
    if (!is(bs, 'BSseq')) {
        stop('bs must be class BSseq. See bsseq::read.bismark() or bsseq::BSseq().')
    }
    if (!is(min_count, 'numeric')) {
        stop('min_count must be an integer.')
    }
    if (!is(max_count, 'numeric')) {
        stop('max_count must be an integer.')
    }
    if (!(min_count < max_count)) {
        stop('min_count not less than max_count')
    }

    cov = bsseq::getCoverage(bs, type = 'Cov')
    meth = bsseq::getCoverage(bs, type = 'M')

    idx = cov < min_count
    cov[idx] = 0
    meth[idx] = 0

    idx = cov > max_count
    cov[idx] = 0
    meth[idx] = 0

    bs = bsseq::BSseq(
        Cov = cov,
        M = meth,
        gr = GenomicRanges::granges(bs),
        pData = pData(bs),
        sampleNames = sampleNames(bs)
    )

    return(bs)
}
