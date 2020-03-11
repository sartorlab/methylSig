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
#' @examples
#' data(BS.cancer.ex, package = 'bsseqData')
#'
#' bs = filter_loci_by_group_coverage(
#'     bs = BS.cancer.ex,
#'     group_column = 'Type',
#'     c('cancer' = 2, 'normal' = 2))
#'
#' small_test = bs[1:50]
#'
#' diff_fit = diff_dss_fit(
#'     bs = small_test,
#'     design = bsseq::pData(bs),
#'     formula = '~ Type')
#'
#' result = diff_dss_test(
#'     diff_fit = diff_fit,
#'     contrast = matrix(c(0,1), ncol = 1)
#' )
#'
#' @export
diff_dss_test = function(
    diff_fit,
    contrast,
    methylation_group_column,
    methylation_groups) {

    #####################################

    # Check missing
    if (missing(diff_fit)) {
        stop('Must pass diff_fit, the result of diff_dss_fit().')
    }
    if (missing(contrast)) {
        stop('Must pass contrast as a matrix.')
    }

    #####################################

    # Check validity of diff_fit
    if (!is(diff_fit, 'list')) {
        stop('diff_fit must be a list.')
    } else {
        if (!all(c('gr', 'design', 'formula', 'X', 'fit') %in% names(diff_fit))) {
            stop('diff_fit must be a list returned from diff_dss_fit() with elements gr, design, formula, X, fit.')
        }
    }

    # Check validity of methylation_group_column
    if (!missing(methylation_group_column)) {
        if (!(is(methylation_group_column, 'character') && length(methylation_group_column) == 1)) {
            stop('methylation_group_column must be a character string.')
        }

        if (!(methylation_group_column %in% colnames(diff_fit$design))) {
            stop(sprintf('methylation_group_column: %s not in column names of diff_fit$design: %s',
                methylation_group_column, paste(colnames(diff_fit$design), collapse = ', ')))
        }
    }

    # Check validity of methylation_groups
    if (!missing(methylation_groups)) {

        if (missing(methylation_group_column)) {
            stop('If methylation_groups is specified, so must methylation_group_column.')
        }

        if (!is(methylation_groups, 'character')) {
            stop('methylation_groups must be a named character vector.')
        }

        if (!all(c('case','control') %in% names(methylation_groups))) {
            stop('methylation_groups must be a named vector with names "case" and "control".')
        }

        if (!all(methylation_groups %in% diff_fit$design[, methylation_group_column])) {
            stop(sprintf('Not all methylation_groups are in methylation_group_column: %s',
                paste(setdiff(methylation_groups, diff_fit$design[, methylation_group_column]), collapse = ', ') ))
        }
    }

    #####################################

    result = DSS::DMLtest.multiFactor(
        DMLfit = diff_fit,
        Contrast = contrast)

    # Create the GRanges return object and harmonize column names with methylSigCalc()
    result_gr = diff_fit$gr
    mcols(result_gr) = result[,c('stat','pvals','fdrs')]
    colnames(mcols(result_gr)) = c('stat','pvalue','fdr')

    return(result_gr)
}
