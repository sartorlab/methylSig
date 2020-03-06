#' Group cytosine / CpG level data into regions based on genomic regions
#'
#' An optional function to aggregate cytosine / CpG level data into regions based on a \code{GRanges} set of genomic regions.
#'
#' @param bs a \code{BSseq} object.
#' @param gr a \code{GRanges} object.
#'
#' @return A \code{BSseq} object with loci of regions matching \code{gr}. Coverage and methylation read count matrices are aggregated by the sums of the cytosines / CpGs in the regions per sample.
#'
#' @examples
#' data(bsseq_stranded, package = 'methylSig')
#' regions = GenomicRanges::GRanges(
#'     seqnames = c('chr1','chr1','chr1'),
#'     ranges = IRanges::IRanges(
#'         start = c(5,35,75),
#'         end = c(30,70,80)
#'     )
#' )
#' tiled = tile_by_regions(bs = bsseq_stranded, gr = regions)
#'
#' @export
tile_by_regions = function(bs, gr) {

    if (missing(bs)) {
        stop('Must pass bs as a BSseq object.')
    }
    if (missing(gr)) {
        stop('Must pass gr as a GRanges object.')
    }
    if (!is(bs, 'BSseq')) {
        stop('bs must be class BSseq.')
    }
    if (!is(gr, 'GRanges')) {
        stop('gr must be class GRanges.')
    }

    #####################################

    cov = bsseq::getCoverage(bs, regions = gr, type = 'Cov', what = 'perRegionTotal')

    if (all(is.na(cov))) {
        stop('No regions overlap between bs and gr')
    }

    meth = bsseq::getCoverage(bs, regions = gr, type = 'M', what = 'perRegionTotal')

    # Set all NA entries to 0 so bsseq::BSseq doesn't throw an error
    # These will likely end up removed in filter_loci_by_group_coverage()
    cov[is.na(cov)] = 0
    meth[is.na(meth)] = 0

    bs = bsseq::BSseq(
        Cov = cov,
        M = meth,
        gr = gr,
        pData = pData(bs),
        sampleNames = sampleNames(bs)
    )

    return(bs)
}
