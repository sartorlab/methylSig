#' Obtain tiled methylation data in non-overlapping continuous windows.
#'
#' This function summarizes methylation data within tiles or user-specified regions. For all CpGs within an intersecting genomic region, the coverage and methylation reads are summed. This is used prior to the \code{\link{methylSigCalc}} function when the user prefers to conduct a tiled analysis instead of a base specific analysis for differential methylation. Tiling may provide higher power to detect significant differences, especially for experiments with low coverage.
#'
#' @param meth A \code{BSseq-class} object, as from \code{methylSigReadData}.
#' @param tiles One of \code{NULL}, a \code{data.frame}, or a \code{GRanges} object. Those CpG sites not belonging any tile will be removed from tiled data.
#' @param win.size An \code{integer} indicating the desired window size in bps. Default is 200. Used only when \code{tiles = NULL}.
#'
#' @return A \code{BSseq-class} object.
#'
#' @examples
#' data(data, package = 'methylSig')
#' methTile = methylSigTile(data)
#' @export
methylSigTile <- function(meth, tiles = NULL, win.size = 200) {
    if(!is(meth, 'BSseq')) {
        stop("'meth' must be a BSseq object.")
    }

    if(mean(width(meth)) != 1) {
        stop("It appears that 'meth' is not CpG resolution. Tiling can only be done on CpG resolution data.")
    }

    # Check for tiles possibilities. Coerce correct seqinfo, and trim().
    if(is.null(tiles)) {
        tiles = tileGenome(seqlengths(meth), tilewidth = 1000, cut.last.tile.in.chrom = TRUE)
        seqinfo(tiles) = seqinfo(meth)
    } else if (is(tiles, 'data.frame')) {
        tiles = makeGRangesFromDataFrame(tiles, seqinfo = seqinfo(meth), keep.extra.columns = FALSE)
        tiles = trim(tiles)
    } else if (is(tiles, 'GRanges')) {
        if(all(is.na(genome(tiles)))) {
            warning("The genome of the GRanges 'tiles' is NA. Coercing to that of 'meth'.")
            seqinfo(tiles) = seqinfo(meth)
        } else if (all(genome(tiles) != genome(meth))) {
            stop("The genome of the GRanges 'tiles' is assigned, but does not match that of 'meth'.")
        }
        tiles = trim(tiles)
        tiles = granges(tiles)
    }

    all_cov = as.matrix(bsseq::getCoverage(meth, type = 'Cov'))
    all_meth = as.matrix(bsseq::getCoverage(meth, type = 'M'))

    overlaps = findOverlaps(meth, tiles, ignore.strand = TRUE)

    if(length(overlaps) > 0) {
        # Prepare for reconstructioon of BSseq object
        tiles = tiles[unique(subjectHits(overlaps))]
        tiled_cov = do.call(rbind, by(all_cov, subjectHits(overlaps), colSums))
        tiled_meth = do.call(rbind, by(all_meth, subjectHits(overlaps), colSums))

        tiled_bsseq = BSseq(M = tiled_meth, Cov = tiled_cov, gr = tiles, pData = pData(meth))

        return(tiled_bsseq)
    } else {
        stop("No regions in 'meth' intersected with regions in 'tiles'.")
    }
}
