#' Group cytosine / CpG level data into regions based on genomic windows
#'
#' An optional function to aggregate cytosine / CpG level data into regions based on a tiling of the genome by \code{win_size}.
#'
#' @param bs a \code{BSseq} object.
#' @param win_size an \code{integer} indicating the size of the tiles. Default is 200bp.
#'
#' @return A \code{BSseq} object with loci consisting of a tiling of the genome by \code{win_size} bp tiles. Coverage and methylation read count matrices are aggregated by the sums of the cytosines / CpGs in the regions per sample.
#'
#' @examples
#' data(bsseq_stranded, package = 'methylSig')
#'
#' tiled = tile_by_windows(bs = bsseq_stranded, win_size = 50)
#'
#' @export
tile_by_windows = function(bs, win_size = 200) {

        if (missing(bs)) {
            stop('Must pass bs as a BSseq object.')
        }
        if (!is(bs, 'BSseq')) {
            stop('bs must be class BSseq.')
        }
        if (!is(win_size, 'numeric')) {
            stop('win_size must be an integer')
        }

        #####################################

        # Determine maximum position per chromosome in use, and add win_size
        seqlevels_in_use = GenomeInfoDb::seqlevelsInUse(bs)
        seqlengths = vapply(seqlevels_in_use, function(chr) {
            gr_tmp = granges(bs)
            chr_length = max(end(gr_tmp[seqnames(gr_tmp) == chr])) + win_size
            return(chr_length)
        }, 1)

        gr = GenomicRanges::tileGenome(
            seqlengths = seqlengths,
            tilewidth = win_size,
            cut.last.tile.in.chrom = TRUE)

        bs = tile_by_regions(bs = bs, gr = gr)

        # To avoid downstream issues with seqinfo mismatches reset lengths to NA
        seqlengths(bs) = NA

        return(bs)
}
