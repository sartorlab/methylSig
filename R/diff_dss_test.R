# Returns contrast as a matrix with one row per column of the design matrix,
# named by coef_names, or stops with the size and names contrast needs
.check_contrast = function(contrast, coef_names) {
    n_coef = length(coef_names)
    contrast_help = sprintf(
        'contrast needs %s rows (or a vector of length %s), one per column of diff_fit$X, in this order:\n%s',
        n_coef, n_coef,
        paste(sprintf('  %s: %s', seq_len(n_coef), coef_names), collapse = '\n'))

    if (!is.numeric(contrast) || length(dim(contrast)) > 2) {
        stop(sprintf('contrast must be a numeric vector or matrix. %s', contrast_help))
    }
    if (is.null(dim(contrast))) {
        contrast = matrix(contrast, ncol = 1, dimnames = list(names(contrast), NULL))
    }
    if (nrow(contrast) != n_coef) {
        stop(sprintf('contrast has %s rows, but %s', nrow(contrast), contrast_help))
    }

    # Named rows are matched to the columns of the design matrix, in any order
    if (is.null(rownames(contrast))) {
        rownames(contrast) = coef_names
    } else if (setequal(rownames(contrast), coef_names) && !anyDuplicated(rownames(contrast))) {
        contrast = contrast[coef_names, , drop = FALSE]
    } else {
        stop(sprintf('The names of contrast (%s) are not the columns of diff_fit$X. %s',
            paste(rownames(contrast), collapse = ', '), contrast_help))
    }

    if (anyNA(contrast)) {
        stop('contrast must not have NA values.')
    }
    if (any(colSums(contrast != 0) == 0)) {
        stop(sprintf('Each column of contrast must have a nonzero value. %s', contrast_help))
    }
    if (qr(contrast)$rank < ncol(contrast)) {
        stop('The columns of contrast must be linearly independent.')
    }

    contrast
}

# Describes the hypothesis a contrast tests, e.g. 'Typenormal = 0', in terms of
# its row names
.describe_contrast = function(contrast) {
    hypotheses = vapply(seq_len(ncol(contrast)), function(j) {
        weights = contrast[, j]
        coefs = rownames(contrast)[weights != 0]
        weights = weights[weights != 0]
        terms = ifelse(abs(weights) == 1, coefs, sprintf('%s * %s', vapply(abs(weights), format, ''), coefs))
        signs = ifelse(weights < 0, '- ', '+ ')
        lhs = paste(signs, terms, sep = '', collapse = ' ')
        sprintf('%s = 0', sub('^\\+ ', '', lhs))
    }, '')
    paste(hypotheses, collapse = ' and ')
}

#' Calculates differential methylation statistics under general experimental design
#'
#' This function is a wrapper for \code{DSS::DMLtest.multiFactor} with the added feature of reporting methylation rates alongside the test results via the \code{methylation_group_column} and \code{methylation_groups} parameters. See documentation below.
#'
#' @param bs a \code{BSseq}, the same used used to create \code{diff_fit}.
#' @param diff_fit a \code{list} object output by \code{diff_dss_fit()}.
#' @param contrast a contrast for hypothesis testing: a numeric vector with one value per column of \code{diff_fit$X}, or a matrix with one row per column of \code{diff_fit$X} and one column per hypothesis. Unnamed values are in the order of the columns of \code{diff_fit$X}, and named values (or row names) are matched to them by name. If the contrast is the wrong size, the error lists the columns of \code{diff_fit$X} in order. A message says what the contrast tests, e.g. \code{Typenormal = 0}.
#' @param methylation_group_column Optionally, a column from \code{diff_fit$design} by which to group samples and capture methylation rates. This column can be a \code{character}, \code{factor}, or \code{numeric}. In the case of \code{numeric} the samples are grouped according to the bottom and top percentiles of the covariate given by \code{covariate_percentiles}, and the mean methylation for each group is calculated. If not a \code{numeric}, use the \code{methylation_groups} parameter to specify case and control.
#' @param methylation_groups Optionally, a named \code{character} vector indicating the \code{case} and \code{control} factors of \code{methylation_group_column} by which to group samples and capture methylation rates. If specified, must also specify \code{methylation_group_column}.
#' @param covariate_percentiles A \code{numeric} vector of two percentiles, from 0 to 100, used when \code{methylation_group_column} is \code{numeric}. Samples with covariate values at or below the first percentile are the case group, and samples at or above the second are the control group. Default \code{c(25, 75)}, the bottom and top 25 percent.
#'
#' @return A \code{GRanges} object containing the following \code{mcols}:
#' \describe{
#'   \item{stat:}{ The test statistic. With a one-column contrast, its sign is the sign of the contrast, e.g. positive when \code{Typenormal > 0}, which need not match the sign of \code{meth_diff}. }
#'   \item{pvalue:}{ The p-value. }
#'   \item{fdr:}{ The Benjamini-Hochberg adjusted p-values using \code{p.adjust(method = 'BH')}. }
#' }
#' If \code{methylation_group_column} is specified, also the following \code{mcols}:
#' \describe{
#'   \item{meth_case:}{ Methylation estimate for case. }
#'   \item{meth_control:}{ Methylation estimate for control. }
#'   \item{meth_diff:}{ The difference \code{meth_case - meth_control}. }
#'   \item{direction:}{ The group for which the locus is hyper-methylated. Note, this is not subject to significance thresholds. }
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
    methylation_groups = NA,
    covariate_percentiles = c(25, 75)) {

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

    # Check validity of covariate_percentiles
    if (!(is(covariate_percentiles, 'numeric') && length(covariate_percentiles) == 2 &&
        !anyNA(covariate_percentiles) && all(covariate_percentiles >= 0 & covariate_percentiles <= 100) &&
        covariate_percentiles[1] < covariate_percentiles[2])) {
        stop('covariate_percentiles must be two increasing numbers from 0 to 100.')
    }

    # Check validity of contrast, and say what it tests, in terms of the
    # columns of diff_fit$X. DSS only checks the number of rows.
    contrast = .check_contrast(contrast, colnames(diff_fit$X))
    message(sprintf('Testing %s', .describe_contrast(contrast)))

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

            # Order of return is covariate_percentiles[1], covariate_percentiles[2]
            # So we want <= quantiles[1] and >= quantiles[2]
            quantiles = quantile(
                x = pdata[, methylation_group_column],
                probs = covariate_percentiles / 100,
                na.rm = TRUE
            )

            case = 'Hyper'
            control = 'Hypo'

            case_idx = which(pdata[, methylation_group_column] <= quantiles[1])
            control_idx = which(pdata[, methylation_group_column] >= quantiles[2])

            # With tied covariate values, both percentiles can be the same value
            if (length(intersect(case_idx, control_idx)) > 0) {
                stop(sprintf('covariate_percentiles %s and %s of methylation_group_column %s put some samples in both groups. Choose percentiles further apart.',
                    covariate_percentiles[1], covariate_percentiles[2], methylation_group_column))
            }

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

        result_gr$direction = ifelse(result_gr$meth_diff >= 0, case, control)

        col_order = c(
            'meth_case',
            'meth_control',
            'meth_diff',
            'direction',
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
