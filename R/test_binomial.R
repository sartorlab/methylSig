#' Differential methylation analysis using binomial model
#'
#' This function calculates differential methylation statistics using a binomial-based approach. See `Warning' message below.
#'
#' This function uses a binomial-based model to calculate differential methylation statistics. It is nearly identical to the \code{methylKit::calculateDiffMeth} function in the \code{methylKit} R package except that only the likelihood ratio test and \code{p.adjust()} with \code{method='BH'} are used to calculate significance levels. It is significantly faster than \code{methylKit::calculateDiffMeth} function.
#'
#' @param meth A \code{BSseq-class} object to calculate differential methylation statistics. See \code{methylSigReadData} for how to read in methylation data.
#' @param comparison The name of the column in \code{pData(meth)} to use for the comparisons, with the correct reference level set.
#' @param min.per.group A vector with two numbers specifying the minimum number of samples required to perform the test for differential methylation. If it is a single number, both groups will use it as the minimum requried number of samples. Default is \code{c(3,3)}.
#'
#' @return \code{GRanges} object containing the differential methylation statistics and locations. \code{p.adjust} with \code{method='BH'} option is used for p-value correction.
#'
#' @section Warning: This function does not take into account the variability among samples in each group being compared.
#'
#' @seealso \code{\link{methylSigCalc}}
#'
#' @examples
#' utils::data(sample_data, package = 'methylSig')
#'
#' myDiff = binomialDiffCalc(meth = meth, comparison = 'DR_vs_DS')
#'
#' @keywords differentialMethylation
#'
#' @export
binomialDiffCalc <- function(
    meth,
    comparison,
    min.per.group=c(3,3)) {

    if(length(min.per.group) == 1) {
        min.per.group = c(min.per.group,min.per.group)
    }

    #####################################
    # Get the group labels, NOTE: THIS ASSUMES CORRECT REFERENCE LEVEL SET
    pdata = bsseq::pData(meth)
    group2 = levels(pdata[, comparison])[2]
    group1 = levels(pdata[, comparison])[1]

    # Determine which rows of pData belong to which group
    # / which columns of Cov and M matrices belong to which group
    group2_idx = which(pdata[,comparison] == group2)
    group1_idx = which(pdata[,comparison] == group1)

    #####################################
    # Determine which sites are valid to test according to min.per.group
    all_cov = as.matrix(bsseq::getCoverage(meth, type = 'Cov'))
    all_meth = as.matrix(bsseq::getCoverage(meth, type = 'M'))

    # Determine which loci satisfy min.per.group
    valid_idx = which(
        base::rowSums(all_cov[, group2_idx] > 0) >= min.per.group[2] & base::rowSums(all_cov[, group1_idx] > 0) >= min.per.group[1]
    )

    #####################################
    # Resize all_cov and all_meth to valid_idx
    all_cov = all_cov[valid_idx,]
    all_meth = all_meth[valid_idx,]

    # Setup required quantities for the logLikRatio calculation
    treads = rowSums(all_cov - all_meth, na.rm = TRUE)
    treads1 = rowSums(all_cov[,group1_idx] - all_meth[,group1_idx], na.rm = TRUE)
    treads2 = rowSums(all_cov[,group2_idx] - all_meth[,group2_idx], na.rm = TRUE)
    creads = rowSums(all_meth, na.rm = TRUE)
    creads1 = rowSums(all_meth[,group1_idx], na.rm = TRUE)
    creads2 = rowSums(all_meth[,group2_idx], na.rm = TRUE)
    cov = rowSums(all_cov, na.rm = TRUE)
    cov1 = rowSums(all_cov[,group1_idx], na.rm=TRUE)
    cov2 = rowSums(all_cov[,group2_idx], na.rm=TRUE)

    logLikRatio = 2 * (creads1 * log(creads1 / cov1 + 1e-100)
                      + treads1 * log(treads1 / cov1 + 1e-100)
                      + creads2 * log(creads2 / cov2 + 1e-100)
                      + treads2 * log(treads2 / cov2 + 1e-100)
                      - creads * log(creads / cov + 1e-100)
                      - treads * log(treads / cov + 1e-100)
                     )

    mu1 = (creads1 / cov1)*100
    mu2 = (creads2 / cov2)*100
    meth.diff = mu2 - mu1

    hyper.direction = ifelse(meth.diff >= 0, group2, group1)

    pvalue = stats::pchisq(logLikRatio, 1, lower.tail=FALSE)
    fdr = stats::p.adjust(pvalue, method = 'BH')

    results = data.frame(
        'logLikRatio' = logLikRatio,
        'mu2' = mu2,
        'mu1' = mu1,
        'meth.diff' = meth.diff,
        'hyper.direction' = hyper.direction,
        'pvalue' = pvalue,
        'fdr' = fdr,
        stringsAsFactors = FALSE
    )



    # Extract the granges of the meth BSseq object and attach the data.frame
    results_gr = GenomicRanges::granges(meth)[valid_idx, ]
    GenomicRanges::mcols(results_gr) = results

    colnames(GenomicRanges::mcols(results_gr)) = c(
        'logLikRatio',
        paste('meth', group2, sep='.'),
        paste('meth', group1, sep='.'),
        'meth.diff',
        'hyper.direction',
        'pvalue',
        'fdr'
    )

    # Add metadata
    results_metadata = list(
        method = 'binomialDiffCalc',
        comparison = comparison,
        min.per.group = min.per.group
    )

    S4Vectors::metadata(results_gr) = results_metadata

    return(results_gr)
}
