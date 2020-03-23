#' Performs model fit for general experimental design
#'
#' This function is a wrapper for \code{DSS::DMLfit.multiFactor}.
#'
#' @param bs a \code{BSseq} object to calculate differential methylation statistics.
#' @param design a \code{data.frame} or \code{DataFrame} for experimental design. Should contain as many rows as there are columns (samples) in \code{bs}, and the order of the rows should match the columns of \code{bs}. If omitted, will default to \code{pData(bs)}.
#' @param formula a formula for the linear model. It should refer to column names from \code{design}. NOTE: The intercept is included by default if omitted. One can omit the intercept with a formula such as \code{'~ 0 + group'}. For clarity, it helps to include the intercept explicitly as in \code{'~ 1 + group'}.
#'
#' @return A \code{list} object with:
#' \describe{
#'     \item{gr:}{ a \code{GRanges} object with loci fit. }
#'     \item{design:}{ the \code{data.frame} input as the experimental design. }
#'     \item{formula:}{ the \code{formula} representing the model. Can be \code{character} or \code{formula}. }
#'     \item{X:}{ the design \code{matrix} used in regression based on the \code{design} and \code{formula}. This should be consulted to determine the appropriate contrast to use in \code{dss_fit_test()}. }
#'     \item{fit:}{ a \code{list} with model fitting results. It has components \code{beta}, the estimated coefficients, and \code{var.beta} the estimated variance/covariance matrix for \code{beta}. }
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
#' @export
diff_dss_fit = function(
    bs,
    design,
    formula) {

    #####################################

    # Check missing
    if (missing(bs)) {
        stop('Must pass bs as a BSseq object.')
    }
    if (missing(formula)) {
        stop('Must pass formula as a character string or formula.')
    }

    #####################################

    # Check types
    if (!is(bs, 'BSseq')) {
        stop('bs must be class BSseq.')
    }
    if (!missing(design)) {
        if (!(any(
                is(design, 'data.frame'),
                is(design, 'DataFrame')
            ))) {
            stop('design must be a data.frame or DataFrame')
        }
    }
    if (!(is(formula, 'character') || is(formula, 'formula'))) {
        stop('formula must be a character string or a formula.')
    }

    #####################################

    if(missing(design)) {
        message('Missing design, defaulting to pData(bs)...')
        design = pData(bs)
    }

    if(is(formula, 'character')) {
        formula = stats::as.formula(formula)
    }

    #####################################

    dss_fit = DSS::DMLfit.multiFactor(
        BSobj = bs,
        design = design,
        formula = formula,
        smoothing = FALSE)

    return(dss_fit)
}
