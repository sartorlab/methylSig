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
#'
#' @return A \code{methylSigData-class} object.
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
#' data = methylSigReadData(
#'     fileList = files,
#'     pData = pData,
#'     assembly = 'hg19',
#'     destranded = TRUE,
#'     maxCount = 500,
#'     minCount = 10,
#'     filterSNPs = TRUE,
#'     num.cores = 4,
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
    fileType = c("cov", "cytosineReport")) {

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
        verbose = TRUE)

    # Filter C > T (+) or G > A (-) SNPs
    # SNPs are 1-based
    # In order to avoid seqinfo() incompatibility issues between the BSseq
    # object and the CT_SNPs_hg19 object
    if(filterSNPs) {
        if(assembly == 'hg19') {
            message('Filtering SNPs')
            data('CT_SNPs_hg19', envir=environment())

            overlaps = findOverlaps(granges(bs), CT_SNPs_hg19, ignore.strand = T)
            snp_invalid_list = queryHits(overlaps)

            bs = bs[-snp_invalid_list]
        } else {
            warning(sprintf('Skipping SNP filtering, %s is not supported.', assembly))
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
    bs = bsseq::BSseq(gr = granges(bs), M = m, Cov = cov, pData = pData, rmZeroCov = TRUE)

    # Add the seqinfo information
    genome_seqinfo = tryCatch({
        GenomeInfoDb::Seqinfo(genome = assembly)
    }, error = function(e){
        # If the genome isn't supported:
        # Get the max within each chrom and add 100bp so the last tile isn't weird
        grl = split(granges(bs), seqnames(bs))
        lengths = sapply(grl, function(gr){
            max(end(gr)) + 100
        })

        g_info = GenomeInfoDb::Seqinfo(
            seqnames = seqnames(seqinfo(bs)),
            seqlengths = lengths,
            isCircular = NA,
            genome = assembly
        )
        return(g_info)
    })
    seqinfo(bs) = genome_seqinfo[seqnames(seqinfo(bs))]

    bs = sort(bs, ignore.strand = TRUE)

    return(bs)
}
