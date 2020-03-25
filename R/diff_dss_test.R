#' Calculates differential methylation statistics under general experimental design
#'
#' This function is a wrapper for \code{DSS::DMLtest.multiFactor} with the added feature of reporting methylation rates alongside the test results via the \code{methylation_group_column} and \code{methylation_groups} parameters. See documentation below.
#'
#' @param bs a \code{BSseq}, the same used used to create \code{diff_fit}.
#' @param diff_fit a \code{list} object output by \code{diff_dss_fit()}.
#' @param contrast a contrast matrix for hypothesis testing. The number of rows should match the number of columns \code{design}. Consult \code{diff_fit$X} to ensure the contrast correponds to the intended test.
#' @param methylation_group_column Optionally, a column from \code{diff_fit$design} by which to group samples and capture methylation rates. This column can be a \code{character}, \code{factor}, or \code{numeric}. In the case of \code{numeric} the samples are grouped according to the top and bottom 25 percentiles of the covariate, and the mean methlyation for each group is calculated. If not a \code{numeric}, use the \code{methylation_groups} parameter to specify case and control.
#' @param methylation_groups Optionally, a named \code{character} vector indicating the \code{case} and \code{control} factors of \code{methylation_group_column} by which to group samples and capture methylation rates. If specified, must also specify \code{methylation_group_column}.
#'
#' @return A \code{GRanges} object containing the following \code{mcols}:
#' \describe{
#'   \item{stat:}{ The test statistic. }
#'   \item{pvalue:}{ The p-value. }
#'   \item{fdr:}{ The Benjamini-Hochberg adjusted p-values using \code{p.adjust(method = 'BH')}. }
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
#'     bs = small_test,
#'     diff_fit = diff_fit,
#'     contrast = matrix(c(0,1), ncol = 1)
#' )
#'
#' result_with_meth = diff_dss_test(
#'     bs = small_test,
#'     diff_fit = diff_fit,
#'     contrast = matrix(c(0,1), ncol = 1),
#'     methylation_group_column = 'Type',
#'     methylation_groups = c('case' = 'cancer', 'control' = 'normal')
#' )
#'
#' @export
diff_dss_test = function(
    bs,
    diff_fit,
    contrast,
    methylation_group_column = NA,
    methylation_groups = NA) {

    #####################################

    # Check missing
    if (missing(bs)) {
        stop('Must pass bs as a BSseq object.')
    }
    if (missing(diff_fit)) {
        stop('Must pass diff_fit, the result of diff_dss_fit().')
    }
    if (missing(contrast)) {
        stop('Must pass contrast as a matrix.')
    }

    #####################################

    # Check validity of bs
    if (!is(bs, 'BSseq')) {
        stop('bs must be class BSseq.')
    }

    # Check validity of diff_fit
    if (!is(diff_fit, 'list')) {
        stop('diff_fit must be a list.')
    } else {
        if (!all(c('gr', 'design', 'formula', 'X', 'fit') %in% names(diff_fit))) {
            stop('diff_fit must be a list returned from diff_dss_fit() with elements gr, design, formula, X, fit.')
        }
    }

    # Check validity of methylation_group_column
    if (!is.na(methylation_group_column)) {
        if (!(is(methylation_group_column, 'character') && length(methylation_group_column) == 1)) {
            stop('methylation_group_column must be a character string.')
        }

        if (!(methylation_group_column %in% colnames(diff_fit$design))) {
            stop(sprintf('methylation_group_column: %s not in column names of diff_fit$design: %s',
                methylation_group_column, paste(colnames(diff_fit$design), collapse = ', ')))
        }
    }

    # Check validity of methylation_groups
    if (!all(is.na(methylation_groups))) {

        if (is.na(methylation_group_column)) {
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

    result_gr = diff_fit$gr
    mcols(result_gr) = result[ ,c('stat','pvals','fdrs')]
    colnames(mcols(result_gr)) = c('stat','pvalue','fdr')

    #####################################

    # If a methylation_group_column is given, retrieve methylation rates
    if (!is.na(methylation_group_column)) {

        # Assign correct case_idx and control_idx based on whether the
        # methylation_group_column is character/factor or numeric
        pdata = diff_fit$design

        if (is(pdata[, methylation_group_column], 'character') || is(pdata[, methylation_group_column], 'factor')) {

            case = methylation_groups['case']
            control = methylation_groups['control']

            case_idx = which(pdata[, methylation_group_column] == case)
            control_idx = which(pdata[, methylation_group_column] == control)

        } else if (is(pdata[, methylation_group_column], 'numeric')) {

            # Order of return is 25%, 75%
            # So we want <= quantiles[1] and >= quantiles[2]
            quantiles = quantile(
                x = pdata[, methylation_group_column],
                probs = c(0.25, 0.75),
                na.rm = TRUE
            )

            case_idx = which(pdata[, methylation_group_column] <= quantiles[1])
            control_idx = which(pdata[, methylation_group_column] >= quantiles[2])

        }

        # Subset bs by what was fit
        result_bs = subsetByOverlaps(bs, diff_fit$gr)

        cov_reads_mat = bsseq::getCoverage(bs, type = 'Cov')
        meth_reads_mat = bsseq::getCoverage(bs, type = 'M')

        # Compute case, control, and methylation difference
        meth_case = (DelayedMatrixStats::rowSums2(
            x = meth_reads_mat,
            cols = case_idx,
            value = TRUE, na.rm = TRUE) / DelayedMatrixStats::rowSums2(
                                                x = cov_reads_mat,
                                                cols = case_idx,
                                                value = TRUE, na.rm = TRUE))

        meth_control = (DelayedMatrixStats::rowSums2(
            x = meth_reads_mat,
            cols = control_idx,
            value = TRUE, na.rm = TRUE) / DelayedMatrixStats::rowSums2(
                                                x = cov_reads_mat,
                                                cols = control_idx,
                                                value = TRUE, na.rm = TRUE))
        meth_diff = meth_case - meth_control

        result_gr$meth_case = round(meth_case * 100, 2)
        result_gr$meth_control = round(meth_control * 100, 2)
        result_gr$meth_diff = round(meth_diff * 100, 2)

        col_order = c(
            'meth_case',
            'meth_control',
            'meth_diff',
            'stat',
            'pvalue',
            'fdr'
        )
        mcols(result_gr) = mcols(result_gr)[, col_order]

    }

    # Remove NA tests and indicate how many failed as in diff_methylsig()
    na_idx = is.na(result_gr$stat)
    if (any(na_idx)) {
        result_gr = result_gr[!na_idx]
        message(sprintf('%s loci were dropped due to insufficient degrees of freedom.', sum(na_idx)))
    }

    return(result_gr)
}
