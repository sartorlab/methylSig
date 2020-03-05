.weight_function <- function(u) (1-u^2)^3

.derivative_phi <- function(phi, local_c, local_t, mu, weight) {
    derivative = 0
    indicator_c = local_c > 0
    indicator_t = local_t > 0
    indicator_ct = local_c + local_t > 0

    if(nrow(local_c) == 1) {
        derivative =
            sum( indicator_c * ( mu * (digamma((mu * phi) + local_c + 1e-100) - digamma(mu * phi + 1e-100)) ) ) +
            sum( indicator_t * ((1 - mu) * (digamma( ((1 - mu) * phi) + local_t + 1e-100) - digamma( ((1-mu) * phi) + 1e-100))) ) -
            sum( indicator_ct * (digamma(phi + local_c + local_t + 1e-100) - digamma(phi)) )
    } else {
        for(g in 1:ncol(local_c)) {
            derivative = derivative +
                sum( indicator_c[,g] * (weight * mu[,g] * (digamma(mu[,g] * phi + local_c[,g] + 1e-100) - digamma(mu[,g] * phi + 1e-100))) ) +
                sum( indicator_t[,g] * (weight * (1 - mu[,g]) * (digamma((1 - mu[,g]) * phi + local_t[,g] + 1e-100) - digamma((1 - mu[,g]) * phi + 1e-100))) ) -
                sum( indicator_ct[,g] * (weight * (digamma(phi + local_c[,g] + local_t[,g] + 1e-100) - digamma(phi))) )
        }
    }

    derivative
}

.derivative_mu <- function(mu, local_c, local_t, phi, weight) {
    derivative = 0
    indicator_c = local_c > 0
    indicator_t = local_t > 0

    if(nrow(local_c) == 1) {
        derivative =
            sum( indicator_c * (digamma(mu * phi + local_c + 1e-100) - digamma(mu * phi + 1e-100)) ) -
            sum( indicator_t * (digamma((1 - mu) * phi + local_t + 1e-100) - digamma((1 - mu) * phi + 1e-100)) )
    } else {
        for(g in 1:ncol(local_c)) {
            derivative = derivative +
                sum( indicator_c[,g] * (weight * (digamma(mu * phi + local_c[,g]+ 1e-100) - digamma(mu * phi + 1e-100))) ) -
                sum( indicator_t[,g] * (weight * (digamma((1 - mu) * phi + local_t[,g] + 1e-100) - digamma((1 - mu) * phi + 1e-100))) )
        }
    }

    derivative
}

.log_likelihood  <- function(mu, phi, local_c, local_t, weight) {
    llik = 0
    indicator_c = local_c > 0
    indicator_t = local_t > 0

    if(nrow(local_c) == 1) {
        llik = llik +
            sum( indicator_c * (lgamma(mu * phi + local_c + 1e-100) - lgamma(mu * phi + 1e-100)) ) +
            sum( indicator_t * (lgamma((1 - mu) * phi + local_t + 1e-100) - lgamma((1 - mu) * phi + 1e-100)) )
    } else {
        for(g in 1:ncol(local_c)) {
            llik = llik +
                sum( indicator_c[,g] * (weight * (lgamma(mu * phi + local_c[,g] + 1e-100) - lgamma(mu * phi + 1e-100))) ) +
                sum( indicator_t[,g] * (weight * (lgamma((1 - mu) * phi + local_t[,g] + 1e-100) - lgamma((1 - mu) + 1e-100))) )
        }
    }

    2*llik
}

#' Calculates differential methylation statistics using a Beta-binomial approach.
#'
#' The function calculates differential methylation statistics between two groups of samples. This is the main function of the methylSig package, and the method most users should use to test for DMCs or DMRs. The function uses a Beta-binomial approach to calculate differential methylation statistics, accounting for variation among samples within each group.
#'
#' The function calculates differential methylation statistics between two groups of samples. The function uses Beta-binomial approach to calculate differential methylation statistics, accounting for variation among samples within each group. Users who wish to tile their data and test for differentially methylated regions (DMRs) instead DMCs should first use the \code{\link{methylSigTile}} function before using this function.
#'
#' @param bs
#' @param group_column
#' @param comparison_groups
#' @param disp_groups
#' @param local_window_size
#' @param local_weight_function A weight kernel function. The input of this function is from -1 to 1. The default is the tri-weight kernel function defined as \code{function(u) = (1-u^2)^3}. Function value and range of parameter for weight function should be from 0 to 1.
#' @param t_approx A \code{logical} value indicating whether to use squared t approximation for the likelihood ratio statistics. Chi-square approximation (\code{t_approx = FALSE}) is recommended when the sample size is large.  Default is \code{TRUE}.
#' @param n_cores An integer denoting how many cores should be used for differential methylation calculations.
#'
#' @return A \code{GRanges} object containing the following \code{mcols}:
#' \describe{
#'   \item{disp_est}{ The disp est. }
#'   \item{log_like_rat}{ The log likelihood ratio. }
#'   \item{df}{ Degrees of freedom used when \code{t_approx = TRUE}. }
#'   \item{meth_est_control}{ Methylation est for control. Groups correspond to the levels in the column used for the group_column in \code{pdata(bs)}. }
#'   \item{meth_est_case}{ Methylation est for case. }
#'   \item{meth_diff}{ The difference \code{meth_est_case - meth_est_control}. }
#'   \item{direction}{ The group for which the CpG/region is hyper-methylated. Groups correspond to the levels in the column used for the group_column in \code{pdata(bs)}. }
#'   \item{pvalue}{ The p-value from the t-test (\code{t_approx = TRUE}) or the Chi-Square test (\code{t_approx = FALSE}). }
#'   \item{fdr}{ The Benjamini-Hochberg adjusted p-values using \code{p.adjust(method = 'BH')}. }
#' }
#'
#' @examples
#' utils::data(sample_data, package = 'methylSig')
#'
#' result = diff_methylsig(
#'     bs = bs,
#'     group_column = 'Type',
#'     comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
#'     disp_groups = c('cancer' = TRUE, 'normal' = TRUE),
#'     local_window_size = 0,
#'     t_approx = TRUE,
#'     n_cores = 1)
#'
#' @export
diff_methylsig = function(bs, group_column, comparison_groups, disp_groups, local_window_size = 0, local_weight_function, t_approx = TRUE, n_cores = 1) {

    # Constants
    min_disp = 1e-6
    min_inverse_disp = 0.001
    max_inverse_disp = max(1/max(min_disp, 1e-6), min_inverse_disp)
    min_meth = 0
    max_meth = 1

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
    if (missing(disp_groups)) {
        stop('Must pass group_column as a logical vector.')
    }

    # Use .weight_function by default, but not in the function definition because
    # this introduces some strange exporting issues. The user really shouldn't
    # have to think about this at all.
    if (missing(local_weight_function)) {
        local_weight_function = .weight_function
    }

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
    if (!is(disp_groups, 'logical')) {
        stop('disp_groups must be a named character vector.')
    }
    if (!is(local_window_size, 'numeric' && length(local_window_size) == 1)) {
        stop('local_window_size must be an integer.')
    }
    if (!is(local_weight_function, 'function')) {
        stop('local_weight_function must be a function.')
    }
    if (!is(t_approx, 'logical')) {
        stop('t_approx must be TRUE/FALSE.')
    }
    if (!is(n_cores, 'numeric' && length(n_cores) == 1)) {
        stop('n_cores must be an integer.')
    }

    # Check valid group_column name
    if (!(group_column %in% colnames(pData(bs)))) {
        stop(sprintf('group_column: %s not in column names of pData(bs): %s',
            group_column, paste(colnames(pData(bs)), collapse = ', ')))
    }

    # Check valid comparison_groups in group_column of pData(bs)
    if (!(all(comparison_groups %in% pData(bs)[, group_column]))) {
        stop(sprintf('Not comparison_groups are in group_column: %s',
            paste(setdiff(comparison_groups, pData(bs)[, group_column]), collapse = ', ') ))
    }

    # Check valid disp_groups in group_column of pData(bs)
    if (!(all(names(disp_groups) %in% pData(bs)[, group_column]))) {
        stop(sprintf('Not all names of disp_groups are in group_column: %s',
            paste(setdiff(names(disp_groups), pData(bs)[, group_column]), collapse = ', ') ))
    }

    #####################################

    # Get the group labels, NOTE: THIS ASSUMES CORRECT REFERENCE LEVEL SET
    pdata = bsseq::pData(bs)
    case = comparison_groups['case']
    control = comparison_groups['control']

    # Determine which rows of pData belong to which group
    # / which columns of Cov and M matrices belong to which group
    case_idx = which(pdata[, group_column] == case)
    control_idx = which(pdata[, group_column] == control)

    # Determine which sample column indexes to use for disp_groups calculation
    if(all(disp_groups)) {
        disp_groups_idx = c(case_idx, control_idx)
    } else if (disp_groups == case) {
        disp_groups_idx = case_idx
    } else if (disp_groups == control) {
        disp_groups_idx = control_idx
    } else {
        stop('"disp_groups" should be one of "both", the name of case, or the name of control')
    }

    #####################################
    num_loci = length(bs)
    gr = granges(bs)

    cov_mat = as.matrix(bsseq::getCoverage(bs, type = 'Cov'))
    meth_mat = as.matrix(bsseq::getCoverage(bs, type = 'M'))

    # est mu per locus within each group. The same value is used for all samples within the same group.
    meth_est = matrix(0, ncol = ncol(bs), nrow = nrow(bs))
    meth_est[, case_idx] = base::rowSums(meth_mat[, case_idx]) / (base::rowSums(cov_mat[, case_idx]) + 1e-100)
    meth_est[, control_idx] = base::rowSums(meth_mat[, control_idx]) / (base::rowSums(cov_mat[, control_idx]) + 1e-100)

    #####################################
    # Go through each valid locus
    result = do.call(rbind, mclapply(seq_along(gr), function(locus_idx){

        if(local_window_size != 0) {
            # NOTE: It is much faster to work with subsets of the result of start()
            # than it is to work with subsets of GRanges.

            # Get the indices which are within the local_window_size, but also limit to 5 CpGs on either side
            local_loci_idx = intersect(
                which(abs(start(gr)[locus_idx] - start(gr)) < local_window_size),
                max(1, locus_idx - 5):min(num_loci, locus_idx + 5))

            if(length(local_loci_idx) == 1) {
                # Do not use local information when there is only one local locus
                local_loci_idx = locus_idx
                local_weights = 1

                # Collect Cov and M matrices for all the loci in the window
                # Rows are loci and columns are samples
                local_cov = matrix(cov_mat[local_loci_idx, ], nrow = 1)
                local_meth = matrix(meth_mat[local_loci_idx, ], nrow = 1)
                local_unmeth = local_cov - local_meth

                # Collect the correct rows of meth_est
                local_meth_est = matrix(meth_est[local_loci_idx, ], nrow = 1)
            } else {
                # We need to scale the loci in the window onto the interval [-1, 1] because
                # that is the domain of the local_weight_function.
                # This is a vector of the distances of the local loci to the loci of interest (domain)
                local_loci_norm = (start(gr)[local_loci_idx] - start(gr)[locus_idx]) / (local_window_size + 1)

                # Calculate the weights
                # Each is a vector of values of the weight function (range)
                local_weights = local_weight_function(local_loci_norm)

                # Collect Cov and M matrices for all the loci in the window
                # Rows are loci and columns are samples
                local_cov = cov_mat[local_loci_idx, ]
                local_meth = meth_mat[local_loci_idx, ]
                local_unmeth = local_cov - local_meth

                # Collect the correct rows of meth_est
                local_meth_est = meth_est[local_loci_idx, ]
            }
        } else {
            # Do not use local information when the local_window_size is 0
            local_loci_idx = locus_idx
            local_weights = 1

            # Collect Cov and M matrices for all the loci in the window
            # Rows are loci and columns are samples
            local_cov = matrix(cov_mat[local_loci_idx, ], nrow = 1)
            local_meth = matrix(meth_mat[local_loci_idx, ], nrow = 1)
            local_unmeth = local_cov - local_meth

            # Collect the correct rows of meth_est
            local_meth_est = matrix(meth_est[local_loci_idx, ], nrow = 1)
        }

        # Compute the degrees of freedom for the locus
        if(disp_groups == 'both') {
            df_subtract = 2
        } else {
            df_subtract = 1
        }
        df = pmax(rowSums(local_cov[, disp_groups_idx, drop = FALSE] > 0) - df_subtract, 0)
        # Compute the degrees of freedom to be used in the test for differential methylation
        df = sum(df * local_weights)

        if(df > 1) {
            ### Common disp_groups calculation
            # This returns a singleton numeric
            if(.derivative_phi(
                phi = max_inverse_disp,
                local_c = local_meth[, disp_groups_idx, drop = FALSE],
                local_t = local_unmeth[, disp_groups_idx, drop = FALSE],
                mu = local_meth_est[, disp_groups_idx, drop = FALSE],
                weight = local_weights) >= 0) {

                disp_est = max_inverse_disp
            } else if(.derivative_phi(
                phi = min_inverse_disp,
                local_c = local_meth[, disp_groups_idx, drop = FALSE],
                local_t = local_unmeth[, disp_groups_idx, drop = FALSE],
                mu = local_meth_est[, disp_groups_idx, drop = FALSE],
                weight = local_weights) <= 0){

                disp_est = min_inverse_disp
            } else {
                disp_est = stats::uniroot(
                    f = .derivative_phi,
                    interval = c(min_inverse_disp, max_inverse_disp),
                    local_meth[, disp_groups_idx, drop = FALSE],
                    local_unmeth[, disp_groups_idx, drop = FALSE],
                    local_meth_est[, disp_groups_idx, drop = FALSE],
                    local_weights)$root
            }

            ### Common group means calculation
            # This returns a numeric vector (control, case, control + case) with the mu est
            group_meth_est_list = list(control_idx, case_idx, c(control_idx, case_idx))
            group_meth_est = rep(0, length(group_meth_est_list))
            for(group_idx in seq_along(group_meth_est_list)) {
                if(sum(local_meth[, group_meth_est_list[[group_idx]], drop = FALSE]) == 0) {
                    # If there are no local C reads, methylation is 0
                    group_meth_est[group_idx] = 0
                } else if (sum(local_unmeth[, group_meth_est_list[[group_idx]], drop = FALSE]) == 0) {
                    # If there are no local T reads, methylation is 1
                    group_meth_est[group_idx] = 1
                } else {
                    # Otherwise, do something fancier
                    group_meth_est[group_idx] = stats::uniroot(
                        f = .derivative_mu,
                        interval = c(min_meth, max_meth),
                        local_meth[, group_meth_est_list[[group_idx]], drop = FALSE],
                        local_unmeth[, group_meth_est_list[[group_idx]], drop = FALSE],
                        disp_est,
                        local_weights)$root
                }
            }

            ### log Likelihood ratio calculation
            log_like_rat =
                .log_likelihood(
                    mu = group_meth_est[1],
                    phi = disp_est,
                    local_c = local_meth[, control_idx, drop = FALSE],
                    local_t = local_unmeth[, control_idx, drop = FALSE],
                    weight = local_weights) +
                .log_likelihood(
                    mu = group_meth_est[2],
                    phi = disp_est,
                    local_c = local_meth[, case_idx, drop = FALSE],
                    local_t = local_unmeth[, case_idx, drop = FALSE],
                    weight = local_weights) -
                .log_likelihood(
                    mu = group_meth_est[3],
                    phi = disp_est,
                    local_c = local_meth[, c(control_idx, case_idx), drop = FALSE],
                    local_t = local_unmeth[, c(control_idx, case_idx), drop = FALSE],
                    weight = local_weights)

            locus_data = c(
                disp_est = disp_est,
                log_like_rat = log_like_rat,
                meth_est_control = group_meth_est[1]*100,
                meth_est_case = group_meth_est[2]*100,
                meth_est_combined = group_meth_est[3]*100,
                df = df + 2)
        } else {
            # Not enough degrees of freedom, return NAs
            locus_data = c(
                disp_est = NA,
                log_like_rat = NA,
                meth_est_control = NA,
                meth_est_case = NA,
                meth_est_combined = NA,
                df = NA)
        }

        return(locus_data)
    }, mc.cores = n_cores))

    # Build GRanges version of result
    result_gr = gr
    S4Vectors::mcols(result_gr) = result

    if(t_approx) {
         result_gr$pvalue = stats::pt(-sqrt(pmax(result_gr$log_like_rat, 0)), result_gr$df) * 2
    } else {
         result_gr$pvalue = stats::pchisq(pmax(result_gr$log_like_rat, 0), 1, lower.tail = F)
    }

    # Set any methylation difference less than 0.01 to 0
    result_gr$meth_diff = (result_gr$meth_est_case - result_gr$meth_est_control)
    result_gr$meth_diff[abs(result_gr$meth_diff) < 0.01] = 0
    result_gr$meth_diff = as.numeric(result_gr$meth_diff)

    result_gr$fdr = stats::p.adjust(result_gr$pvalue, method = 'BH')

    result_gr$direction = ifelse(result_gr$meth_diff >= 0, case, control)

    col_order = c(
        'disp_est',
        'log_like_rat',
        'df',
        'meth_est_case',
        'meth_est_control',
        'meth_diff',
        'direction',
        'pvalue',
        'fdr'
    )

    S4Vectors::mcols(result_gr) = S4Vectors::mcols(result_gr)[, col_order]

    results_metadata = list(
        method = 'diff_methylsig',
        group_column = group_column,
        disp_groups = disp_groups,
        local_window_size = local_window_size,
        local_weight_function,
        t_approx = t_approx
    )

    S4Vectors::metadata(result_gr) = results_metadata

    ############################################################################

    return(result_gr)
}
