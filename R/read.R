#' Read methylation score files to make a 'BSseq' object.
#'
#' This function reads files created by the Bismark Methylation Extractor, and outputs a \code{BSseq} object.
#'
#' @param fileList Files to be read. These can be \code{cov} or \code{cytosine_reports} from the Bismark Methylation Extractor. See \code{fileType} for details.
#' @param pData A \code{data.frame} containing phenotype information for the samples in \code{fileList}. The \code{row.names} attribute of the \code{data.frame} should match the \code{Sample_Names}. See example below.
#' @param assembly The genome assembly used for alignment. e.g. \code{hg19}, \code{mm10}, etc.
#' @param destranded A logical value indicating whether to destrand the reverse to forward strand. If TRUE, the reads from both will be combined. Default is TRUE.
#' @param maxCount A number indicating the maximum coverage count to be included.
#' @param minCount A number indicating the minimum coverage count to be included.
#' @param filterSNPs A logical value indicating whether or not to filter out C > T SNPs based on the 1000 Genomes Project. NOTE: Only supported when \code{assembly = 'hg19'}.
#' @param num.cores Number of cores to be used in reading files. Default is 1.
#' @param fileType The format of the input file. Either \code{cov} or \code{cytosineReport}. One of the outputs of the Bismark Methylation Extractor.
#' @param verbose A logical value indicating whether \code{bsseq::read.bismark} shoud print progress. Default TRUE.
#'
#' @return A \code{BSseq-class} object.
#'
#' @seealso \code{\link{methylSigCalc}}
#'
#' @examples
#' files = c(
#'     system.file('extdata', 'MDAMB_231_1DR.txt.gz', package='methylSig'),
#'     system.file('extdata', 'MDAMB_231_1DS.txt.gz', package='methylSig'),
#'     system.file('extdata', 'MDAMB_231_2DR.txt.gz', package='methylSig'),
#'     system.file('extdata', 'MDAMB_231_2DS.txt.gz', package='methylSig'),
#'     system.file('extdata', 'MDAMB_231_3DR.txt.gz', package='methylSig'),
#'     system.file('extdata', 'MDAMB_231_3DS.txt.gz', package='methylSig'))
#'
#' sample.ids = basename(files)
#' sample.ids = gsub('.txt.gz', '', sample.ids)
#'
#' pData = data.frame(
#'     Sample_Names = sample.ids,
#'     DR_vs_DS = relevel(factor(c('DR','DS','DR','DS','DR','DS')), ref = 'DS'),
#'     row.names = sample.ids,
#'     stringsAsFactors = FALSE)
#'
#' meth = methylSigReadData(
#'     fileList = files,
#'     pData = pData,
#'     assembly = 'hg19',
#'     destranded = TRUE,
#'     maxCount = 500,
#'     minCount = 10,
#'     filterSNPs = TRUE,
#'     num.cores = 1,
#'     fileType = 'cytosineReport')
#'
#' @export
methylSigReadData = function(
    fileList,
    pData,
    assembly = NA,
    destranded = TRUE,
    maxCount = 500, minCount = 10,
    filterSNPs = FALSE,
    num.cores = 1,
    fileType = c("cov", "cytosineReport"),
    verbose = TRUE) {

    # NOTE: The cytosine report is 1-based, and GRanges is also 1-based. The result,
    # of bsseq::read.bismark() is 1-based. We are 1-based y'all!

    fileType = match.arg(fileType)

    # Read
    bs = bsseq::read.bismark(
        files = fileList,
        sampleNames = rownames(pData),
        rmZeroCov = TRUE,
        strandCollapse = destranded,
        fileType = fileType,
        mc.cores = num.cores,
        verbose = verbose)

    # Assign Seqinfo to bs
    # If the assembly is supported, then downstream functions such as tiling and
    # annotation should work without issue, assuming that the annotations for
    # the same genome don't somehow have different seqlengths. This shouldn't
    # happen on-the-fly because GenomeInfoDb::Seqinfo is used. BUT, if saved
    # annotations are used that are old, GenomeInfoDb::Seqinfo could result in
    # mismatched lengths between Bioc versions

    # If the assembly is not supported, we have to be careful about how the
    # Seqinfo is defined. If we try to define it as the end of each present
    # chromosome + some width to allow for tiling, we run the risk of not being
    # able to perform custom annotations if the annotation Seqinfo doesn't match.

    # If the Seqinfo is left blank, we cannot tile downstream, warn the user.
    if(!is.na(assembly)) {
        genome_seqinfo = tryCatch({
            GenomeInfoDb::Seqinfo(genome = assembly)
        }, error = function(e) {
            warning(sprintf('The assembly %s is not supported by GenomeInfoDb::fetchExtendedChromInfoFromUCSC, and the resulting BSseq object from this function will have the trivial seqinfo attributes. In order to use the methylSigTile function, you must create a custom GenomeInfoDb::Seqinfo and assign it to the result of this function. If you would like to use annotation functions downstream, the seqinfo for custom annotations MUST BE THE SAME as what is assigned to the result this function.', assembly))
            seqinfo(bs)
        })
    } else {
        warning('Leaving assembly as NA will give the resulting BSseq object from this function the trivial seqinfo attributes. In order to use the methylSigTile function downstream, you MUST create a custom GenomeInfoDb::Seqinfo and assign it to the result of this function. If you would like to use annotation functions downstream, the seqinfo for custom annotations MUST BE THE SAME as what is assigned to the result of this function.')

        genome_seqinfo = seqinfo(bs)
    }
    seqinfo(bs) = merge(seqinfo(bs), genome_seqinfo)

    # Filter C > T (+) or G > A (-) SNPs
    # SNPs are 1-based
    # In order to avoid seqinfo() incompatibility issues between the BSseq
    # object and the CT_SNPs_hg19 object
    if(filterSNPs) {
        if(!is.na(assembly) && assembly == 'hg19') {
            message('Filtering SNPs')
            utils::data('CT_SNPs_hg19', envir=environment())
            CT_SNPs_hg19 = get('CT_SNPs_hg19')

            overlaps = GenomicRanges::findOverlaps(GenomicRanges::granges(bs), CT_SNPs_hg19, ignore.strand = T)
            snp_invalid_list = S4Vectors::queryHits(overlaps)

            bs = bs[-snp_invalid_list]
        } else {
            message(sprintf('Skipping SNP filtering, genome %s is not supported.', assembly))
        }
    }

    # Filter for maxCount and minCount
    # Must change both Cov and M because BSseq has some sanity checks
    # where 0 <= M <= Cov
    cov  = as.matrix(bsseq::getCoverage(BSseq = bs, type = 'Cov'))
    m = as.matrix(bsseq::getCoverage(BSseq = bs, type = 'M'))
    for(j in 1:ncol(cov)) {
        count_idx = which(cov[,j] < minCount | cov[,j] > maxCount)
        cov[count_idx, j] = 0
        m[count_idx, j] = 0
    }

    # Rebuild the BSseq object after altering the Cov and M counts
    bs = bsseq::BSseq(gr = GenomicRanges::granges(bs), M = m, Cov = cov, pData = pData, rmZeroCov = TRUE)

    bs = sort(bs, ignore.strand = TRUE)

    bs_metadata = list(
        files = fileList,
        assembly = assembly,
        destranded = destranded,
        maxCount = maxCount,
        minCount = minCount,
        filterSNPs = filterSNPs,
        fileType = fileType,
        tile = FALSE,
        tiles = NULL,
        win.size = NULL
    )

    S4Vectors::metadata(bs) = bs_metadata

    return(bs)
}
