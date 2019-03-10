#' Obtain tiled methylation data in non-overlapping continuous windows.
#'
#' This function summarizes methylation data within tiles or user-specified regions. For all CpGs within an intersecting genomic region, the coverage and methylation reads are summed. This is used prior to the \code{\link{methylSigCalc}} function when the user prefers to conduct a tiled analysis instead of a base specific analysis for differential methylation. Tiling may provide higher power to detect significant differences, especially for experiments with low coverage.
#'
#' @param meth A \code{BSseq-class} object, as from \code{methylSigReadData}.
#' @param tiles One of \code{NULL}, a \code{data.frame}, or a \code{GRanges} object. If not \code{NULL}, the regions should be non-overlapping. Those CpG sites not belonging to any tile will be removed from tiled data.
#' @param win.size An \code{integer} indicating the desired window size in bps. Default is 200. Used only when \code{tiles = NULL}.
#'
#' @return A \code{BSseq-class} object.
#'
#' @examples
#' utils::data(sample_data, package = 'methylSig')
#' methTile = methylSigTile(meth, tiles = NULL, win.size = 200)
#'
#' @export
methylSigTile <- function(meth, tiles = NULL, win.size = 200) {
    if(!is(meth, 'BSseq')) {
        stop("'meth' must be a BSseq object.")
    }

    if(mean(BiocGenerics::width(meth)) != 1) {
        stop("It appears that 'meth' is not CpG resolution. Tiling can only be done on CpG resolution data.")
    }

    # Check for tiles possibilities
    if(is.null(tiles)) {
        # If the seqlengths aren't defined, remind the user to create a custom GenomeInfoDb::Seqinfo and assign it to meth
        # if(any(is.na(GenomeInfoDb::seqlengths(meth)))) {
        #     stop("The seqinfo for 'meth' is ill-defined, with seqlengths being NA. In order to use the methylSigTile function, you should create a custom GenomeInfoDb::Seqinfo and assign it to 'meth'.")
        # }

        seqlevels_in_use = seqlengths(meth)[seqlevelsInUse(meth)]
        tiles = GenomicRanges::tileGenome(seqlevels_in_use, tilewidth = win.size, cut.last.tile.in.chrom = TRUE)
        seqinfo(tiles) = merge(seqinfo(tiles), seqinfo(meth))
    } else if (is(tiles, 'data.frame')) {
        tiles = GenomicRanges::makeGRangesFromDataFrame(tiles, keep.extra.columns = FALSE)
        seqinfo(tiles) = merge(seqinfo(tiles), seqinfo(meth))
    } else if (is(tiles, 'GRanges')) {
        tiles = GenomicRanges::granges(tiles)
        seqinfo(tiles) = merge(seqinfo(tiles), seqinfo(meth))
    }

    # Subset tiles based on findOverlaps to save some work downstream
    overlaps = GenomicRanges::findOverlaps(query = tiles, subject = meth)
    tile_idx = S4Vectors::queryHits(overlaps)
    tiles = tiles[unique(tile_idx)]

	tiled_M = as.matrix(bsseq::getCoverage(BSseq = meth, regions = tiles, what = "perRegionTotal", type = 'M'))
	tiled_M[is.na(tiled_M)] = 0
	tiled_Cov = as.matrix(bsseq::getCoverage(BSseq = meth, regions = tiles, what = "perRegionTotal", type = 'Cov'))
	tiled_Cov[is.na(tiled_Cov)] = 0

	tiled_bsseq = bsseq::BSseq(
        gr = tiles,
        M = tiled_M,
        Cov = tiled_Cov,
        pData = bsseq::pData(meth),
        sampleNames = rownames(bsseq::pData(meth)),
        rmZeroCov = TRUE)

    S4Vectors::metadata(tiled_bsseq) = S4Vectors::metadata(meth)
    S4Vectors::metadata(tiled_bsseq)$tile = TRUE
    S4Vectors::metadata(tiled_bsseq)$tiles = ifelse(is.null(tiles), 'windows', 'custom')
    S4Vectors::metadata(tiled_bsseq)$win.size = win.size
    S4Vectors::metadata(tiled_bsseq)$cpgs.per.tile = as.numeric(table(tile_idx))

    return(tiled_bsseq)
}
