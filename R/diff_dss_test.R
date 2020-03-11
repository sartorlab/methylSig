#' Calculates differential methylation statistics under general experimental design
#'
#' This function is a wrapper for \code{DSS::DMLfit.multiFactor} and \code{DSS::DMLtest.multiFactor} with the added feature of reporting methylation rates alongside the test results via the \code{methylation_group_column} and \code{methylation_groups} parameters. See documentation below.
#'
#' @param diff_fit a \code{BSseq-class} object to calculate differential methylation statistics.
#' @param contrast a contrast matrix for hypothesis testing. The number of rows should match the number of columns \code{design}.
#' @param methylation_group_column a thing.
#' @param methylation_groups a thing.
#'
#' @return A \code{GRanges} object containing the following \code{mcols}:
#' \describe{
#'   \item{stat}{ The test statistic }
#'   \item{pvalue}{ p-values }
#'   \item{fdr}{ Degrees of freedom used when \code{T.approx = TRUE}. }
#' }
#'
#' @export
diff_dss_test = function(
    diff_fit,
    contrast,
    methylation_group_column,
    methylation_groups) {

    result = DSS::DMLtest.multiFactor(
        DMLfit = diff_fit,
        Contrast = contrast)

    # Create the GRanges return object and harmonize column names with methylSigCalc()
    result_gr = diff_fit$gr
    mcols(result_gr) = result[,c('stat','pvals','fdrs')]
    colnames(mcols(result_gr)) = c('stat','pvalue','fdr')

    return(result_gr)
}
