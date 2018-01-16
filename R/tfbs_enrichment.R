#' Perform transcription factor enrichment test among differentially methylated cytosines or regions
#'
#' This function tests for enriched transcription binding sites among differentially methylated sites or regions using a binomial test.
#'
#' Likelihood ratio test is used based on the binomial distribution.
#'
#' @param myDiff \code{GRanges} object resulting from \code{methylSigCalc} that contains all CpG sites that are tested for differential methylation.
#' @param dmcList A \code{logical} of the same length as \code{myDiff} defining the DMCs or DMRs.
#' @param tfbsInfo A \code{GRanges} object of the genomic regions representing peaks. The \code{name} column should indicate which TF the peak is for.
#'
#' @return A \code{data.frame} whose \code{rownames} are inherited from the \code{name} column of the input BED file, and whose columns are:
#' \describe{
#'   \item{n_total_by_tf}{ The number of tested CpGs in a TFBS for a TF. }
#'   \item{n_dmc_by_tf}{ The number of DM CpGs in a TFBS for a TF. }
#'   \item{N_total}{ The total number of tested CpGs in a TFBS across all the TFs. }
#'   \item{N_dmc}{ The total number of DM CpGs in a TFBS across all the TFs. }
#'   \item{p_total}{ \code{n_total_by_tf} / \code{N_total}, used in the likelihood calculation. }
#'   \item{p_dmc}{ \code{n_dmc_by_tf} / \code{N_dmc}, used in the likelihood calculation. }
#'   \item{logLik}{ The log-likelihood based on the binomial distribution. }
#'   \item{pvalue}{ The p-value from the likelihood ratio test. }
#' }
#'
#' @examples
#' data(data, package = 'methylSig')
#'
#' dmcList = msig_cpgs$fdr < 0.05 & abs(msig_cpgs$meth.diff) > 25
#'
#' methylSig.tfbsEnrichTest(myDiff = msig_cpgs, dmcList = dmcList, tfbsInfo = tfbs)
#'
#' @export
methylSig.tfbsEnrichTest <- function(myDiff, dmcList, tfbsInfo) {
    # NOTE: All notation is relative to the paper

    # Create GRangesList of tfbsInfo based on 'name'
    tfbs_by_tf = split(tfbs, f = tfbs$name)

    # Subset myDiff by the dmcList
    dmcs = myDiff[dmcList]

    # Overlap tested CpGs and DMCs with TFBSs, respectively
    tested_overlaps = GenomicRanges::findOverlaps(tfbs, myDiff)
    dmc_overlaps = GenomicRanges::findOverlaps(tfbs, dmcs)

    # Determine the total number of CpGs and DMCs in TFBSs, respectively
    N_total = length(unique(S4Vectors::subjectHits(tested_overlaps)))
    N_dmc = length(unique(S4Vectors::subjectHits(dmc_overlaps)))

    # Per TF overlaps
    by_tf = do.call(rbind, lapply(tfbs_by_tf, function(tf) {
        tested_overlaps_tf = findOverlaps(tf, myDiff)
        dmc_overlaps_tf = findOverlaps(tf, dmcs)

        if(length(tested_overlaps_tf) == 0) {
            return(c('n_total_by_tf' = 0, 'n_dmc_by_tf' = 0))
        } else {
            return(c(
                'n_total_by_tf' = length(unique(S4Vectors::subjectHits(tested_overlaps_tf))),
                'n_dmc_by_tf' = length(unique(S4Vectors::subjectHits(dmc_overlaps_tf)))
            ))
        }
    }))
    by_tf = data.frame(by_tf)

    n_total_by_tf = by_tf$n_total_by_tf
    n_dmc_by_tf = by_tf$n_dmc_by_tf

    p_total = n_total_by_tf / N_total
    p_dmc = n_dmc_by_tf / N_dmc

    # Spacing is to help figure out what's grouped together
    logLik = 2    *    ( n_dmc_by_tf * log( pmax(p_dmc, 1e-100)  /  p_total )   +   (N_dmc - n_dmc_by_tf) * log( pmax(1 - p_dmc, 1e-100)  /  (1 - p_total) ) )

    pvalue = pchisq(logLik, 1, lower.tail=FALSE)

    by_tf$N_total = N_total
    by_tf$N_dmc = N_dmc
    by_tf$p_total = p_total
    by_tf$p_dmc = p_dmc
    by_tf$logLik = logLik
    by_tf$pvalue = pvalue

    # Why this step?
    by_tf[p_dmc < p_total, 'pvalue'] = 1

    by_tf = subset(by_tf, !is.na(by_tf$pvalue))

    return(by_tf)
}
