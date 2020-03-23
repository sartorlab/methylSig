#' Differential methylation analysis using binomial model
#'
#' This function calculates differential methylation statistics using a binomial-based approach. See `Warning' message below.
#'
#' This function uses a binomial-based model to calculate differential methylation statistics. It is nearly identical to the \code{methylKit::calculateDiffMeth} function in the \code{methylKit} R package except that only the likelihood ratio test and \code{p.adjust(..., method='BH')} are used to calculate significance levels. It is significantly faster than \code{methylKit::calculateDiffMeth} function.
#'
#' @param bs A \code{BSseq-class} object to calculate differential methylation statistics. See \code{methylSigReadData} for how to read in methylation data.
#' @param group_column a \code{character} string indicating the column of \code{pData(bs)} to use for determining group membership.
#' @param comparison_groups a named \code{character} vector indicating the \code{case} and \code{control} factors of \code{group_column} for the comparison.
#'
#' @return A \code{GRanges} object containing the following \code{mcols}:
#' \describe{
#'   \item{meth_case:}{ Methylation estimate for case. }
#'   \item{meth_control:}{ Methylation estimate for control. }
#'   \item{meth_diff:}{ The difference \code{meth_case - meth_control}. }
#'   \item{direction:}{ The group for which the lcous is hyper-methylated. Note, this is not subject to significance thresholds. }
#'   \item{pvalue:}{ The p-value from the t-test (\code{t_approx = TRUE}) or the Chi-Square test (\code{t_approx = FALSE}). }
#'   \item{fdr:}{ The Benjamini-Hochberg adjusted p-values using \code{p.adjust(method = 'BH')}. }
#'   \item{log_lik_ratio:}{ The log likelihood ratio. }
#' }
#'
#' @section Warning: This function does not take into account the variability among samples in each group being compared.
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
#' diff_gr = diff_binomial(
#'     bs = small_test,
#'     group_column = 'Type',
#'     comparison_groups = c('case' = 'cancer', 'control' = 'normal'))
#'
#' @export
diff_binomial = function(
    bs,
    group_column,
    comparison_groups) {

    #####################################

    # Check missing
    if (missing(bs)) {
        stop('Must pass bs as a BSseq object.')
    }
    if (missing(group_column)) {
        stop('Must pass group_column as a character string.')
    }
    if (missing(comparison_groups)) {
        stop('Must pass comparison_groups as a named character vector with names "case" and "control".')
    }

    #####################################

    # Check types
    if (!is(bs, 'BSseq')) {
        stop('bs must be class BSseq.')
    }
    if (!(is(group_column, 'character') && length(group_column) == 1)) {
        stop('group_column must be a character string.')
    }
    if (!is(comparison_groups, 'character')) {
        stop('comparison_groups must be a named character vector.')
    }

    #####################################

    # Check valid group_column name
    if (!(group_column %in% colnames(pData(bs)))) {
        stop(sprintf('group_column: %s not in column names of pData(bs): %s',
            group_column, paste(colnames(pData(bs)), collapse = ', ')))
    }

    # Check valid comparison_groups values in group_column of pData(bs)
    if (!all(comparison_groups %in% pData(bs)[, group_column])) {
        stop(sprintf('Not all comparison_groups are in group_column: %s',
            paste(setdiff(comparison_groups, pData(bs)[, group_column]), collapse = ', ') ))
    }

    # Check valid comparison_groups names
    if (!all(c('case','control') %in% names(comparison_groups))) {
        stop('comparison_groups vector must be a named vector with names "case" and "control".')
    }

    #####################################

    case = comparison_groups['case']
    control = comparison_groups['control']

    # Rows of pdata and columns of bs
    pdata = bsseq::pData(bs)
    case_idx = which(pdata[, group_column] == case)
    control_idx = which(pdata[, group_column] == control)

    #####################################

    gr = granges(bs)

    cov_mat = as.matrix(bsseq::getCoverage(bs, type = 'Cov'))
    meth_mat = as.matrix(bsseq::getCoverage(bs, type = 'M'))

    # Determine which sites are valid to test according to min.per.group
    cov_mat = as.matrix(bsseq::getCoverage(bs, type = 'Cov'))
    meth_mat = as.matrix(bsseq::getCoverage(bs, type = 'M'))

    # Setup required quantities for the log_lik_ratio calculation
    unmeth_reads = rowSums(cov_mat - meth_mat, na.rm = TRUE)
    unmeth_reads_control = rowSums(cov_mat[,control_idx] - meth_mat[,control_idx], na.rm = TRUE)
    unmeth_reads_case = rowSums(cov_mat[,case_idx] - meth_mat[,case_idx], na.rm = TRUE)

    meth_reads = rowSums(meth_mat, na.rm = TRUE)
    meth_reads_control = rowSums(meth_mat[,control_idx], na.rm = TRUE)
    meth_reads_case = rowSums(meth_mat[,case_idx], na.rm = TRUE)

    cov = rowSums(cov_mat, na.rm = TRUE)
    cov_control = rowSums(cov_mat[,control_idx], na.rm=TRUE)
    cov_case = rowSums(cov_mat[,case_idx], na.rm=TRUE)

    log_lik_ratio = 2 * (meth_reads_control * log(meth_reads_control / cov_control + 1e-100)
                      + unmeth_reads_control * log(unmeth_reads_control / cov_control + 1e-100)
                      + meth_reads_case * log(meth_reads_case / cov_case + 1e-100)
                      + unmeth_reads_case * log(unmeth_reads_case / cov_case + 1e-100)
                      - meth_reads * log(meth_reads / cov + 1e-100)
                      - unmeth_reads * log(unmeth_reads / cov + 1e-100)
                     )

    meth_control = round((meth_reads_control / cov_control) * 100, 2)
    meth_case = round((meth_reads_case / cov_case) * 100, 2)
    meth_diff = round(meth_case - meth_control, 2)

    direction = ifelse(meth_diff >= 0, case, control)

    pvalue = stats::pchisq(log_lik_ratio, 1, lower.tail=FALSE)
    fdr = stats::p.adjust(pvalue, method = 'BH')

    results = data.frame(
        'meth_case' = meth_case,
        'meth_control' = meth_control,
        'meth_diff' = meth_diff,
        'direction' = direction,
        'pvalue' = pvalue,
        'fdr' = fdr,
        'log_lik_ratio' = log_lik_ratio,
        stringsAsFactors = FALSE
    )

    result_gr = granges(bs)
    mcols(result_gr) = results

    return(result_gr)
}
