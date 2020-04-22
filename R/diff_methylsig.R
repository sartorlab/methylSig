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
        for(g in seq(ncol(local_c))) {
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
        for(g in seq(ncol(local_c))) {
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
        for(g in seq(ncol(local_c))) {
            llik = llik +
                sum( indicator_c[,g] * (weight * (lgamma(mu * phi + local_c[,g] + 1e-100) - lgamma(mu * phi + 1e-100))) ) +
                sum( indicator_t[,g] * (weight * (lgamma((1 - mu) * phi + local_t[,g] + 1e-100) - lgamma((1 - mu) + 1e-100))) )
        }
    }

    2*llik
}

#' Calculates differential methylation statistics using a Beta-binomial approach
#'
#' The function calculates differential methylation statistics between two groups of samples using a beta-binomial approach to calculate differential methylation statistics, accounting for variation among samples within each group. The function can be applied to a \code{BSseq} object subjected to \code{filter_loci_by_coverage()}, \code{filter_loci_by_snps()}, \code{filter_loci_by_group_coverage()} or any combination thereof. Moreover, the function can be applied to a \code{BSseq} object which has been tiled with \code{tile_by_regions()} or \code{tile_by_windows()}.
#'
#' @param bs a \code{BSseq} object.
#' @param group_column a \code{character} string indicating the column of \code{pData(bs)} to use for determining group membership.
#' @param comparison_groups a named \code{character} vector indicating the \code{case} and \code{control} factors of \code{group_column} for the comparison.
#' @param disp_groups a named \code{logical} vector indicating the whether to use \code{case}, \code{control}, or both to estimate the dispersion.
#' @param local_window_size an \code{integer} indicating the size of the window for use in determining local information to improve mean and dispersion parameter estimations. In addition to a the distance constraint, a maximum of 5 loci upstream and downstream of the locus are used. The default is \code{0}, indicating no local information is used.
#' @param local_weight_function a weight kernel function. The default is the tri-weight kernel function defined as \code{function(u) = (1-u^2)^3}. The domain of any given weight function should be [-1,1], and the range should be [0,1].
#' @param t_approx a \code{logical} value indicating whether to use squared t approximation for the likelihood ratio statistics. Chi-square approximation (\code{t_approx = FALSE}) is recommended when the sample size is large.  Default is \code{TRUE}.
#' @param n_cores an \code{integer} denoting how many cores should be used for differential methylation calculations.
#'
#' @return A \code{GRanges} object containing the following \code{mcols}:
#' \describe{
#'   \item{meth_case:}{ Methylation estimate for case. }
#'   \item{meth_control:}{ Methylation estimate for control. }
#'   \item{meth_diff:}{ The difference \code{meth_case - meth_control}. }
#'   \item{direction:}{ The group for which the locus is hyper-methylated. Note, this is not subject to significance thresholds. }
#'   \item{pvalue:}{ The p-value from the t-test (\code{t_approx = TRUE}) or the Chi-Square test (\code{t_approx = FALSE}). }
#'   \item{fdr:}{ The Benjamini-Hochberg adjusted p-values using \code{p.adjust(method = 'BH')}. }
#'   \item{disp_est:}{ The dispersion estimate. }
#'   \item{log_lik_ratio:}{ The log likelihood ratio. }
#'   \item{df:}{ Degrees of freedom used when \code{t_approx = TRUE}. }
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
#' small_test = bs[seq(50)]
#'
#' diff_gr = diff_methylsig(
#'     bs = small_test,
#'     group_column = 'Type',
#'     comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
#'     disp_groups = c('case' = TRUE, 'control' = TRUE),
#'     local_window_size = 0,
#'     t_approx = TRUE,
#'     n_cores = 1)
#'
#' @export
diff_methylsig = function(
    bs,
    group_column,
    comparison_groups,
    disp_groups,
    local_window_size = 0,
    local_weight_function,
    t_approx = TRUE,
    n_cores = 1) {

    # Constants
    min_disp = 1e-6
    min_inverse_disp = 0.001
    max_inverse_disp = max(1/max(min_disp, 1e-6), min_inverse_disp)
    min_meth = 0
    max_meth = 1

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
    if (missing(disp_groups)) {
        stop('Must pass disp_groups as a logical vector.')
    }

    # Use .weight_function by default, but not in the function definition because
    # this introduces some strange exporting issues. The user really shouldn't
    # have to think about this at all.
    if (missing(local_weight_function)) {
        local_weight_function = .weight_function
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
    if (!is(disp_groups, 'logical')) {
        stop('disp_groups must be a named logical vector.')
    }
    if (!(is(local_window_size, 'numeric') && length(local_window_size) == 1)) {
        stop('local_window_size must be an integer.')
    }
    if (!is(local_weight_function, 'function')) {
        stop('local_weight_function must be a function.')
    }
    if (!is(t_approx, 'logical')) {
        stop('t_approx must be TRUE/FALSE.')
    }
    if (!(is(n_cores, 'numeric') && length(n_cores) == 1)) {
        stop('n_cores must be an integer.')
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

    # Check valid disp_groups names
    if (!all(c('case','control') %in% names(disp_groups))) {
        stop('disp_groups vector must be a named vector with names "case" and "control".')
    }

    # Check valid disp_groups values
    if (!any(disp_groups)) {
        stop('disp_groups must be a named logical vector with at least one TRUE value corresponding to group name for case or control.')
    }

    # Check for invalid local_window_size == 0 && regions state
    # Cannot use local information on region-resolution data, that's the point of tiling
    if (local_window_size > 0 && median(width(bs)) > 2) {
        stop('Cannot use local information on region-resolution data. Detected local_window_size > 0 and median width of loci > 2')
    }

    #####################################

    case = comparison_groups['case']
    control = comparison_groups['control']

    # Rows of pdata and columns of bs
    pdata = bsseq::pData(bs)
    case_idx = which(pdata[, group_column] == case)
    control_idx = which(pdata[, group_column] == control)

    if(all(disp_groups)) {
        disp_groups_idx = c(case_idx, control_idx)
    } else if (disp_groups['case'] & !disp_groups['control']) {
        disp_groups_idx = case_idx
    } else if (!disp_groups['case'] & disp_groups['control']) {
        disp_groups_idx = control_idx
    }

    #####################################

    num_loci = length(bs)
    gr = granges(bs)

    cov_mat = as.matrix(bsseq::getCoverage(bs, type = 'Cov'))
    meth_mat = as.matrix(bsseq::getCoverage(bs, type = 'M'))

    # Estimate meth per locus within each group. The same value is used for all samples within the same group.
    # Note, the approach is to sum reads over all samples per group per locus
    meth_est = matrix(0, ncol = ncol(bs), nrow = nrow(bs))
    meth_est[, case_idx] = base::rowSums(meth_mat[, case_idx]) / (base::rowSums(cov_mat[, case_idx]) + 1e-100)
    meth_est[, control_idx] = base::rowSums(meth_mat[, control_idx]) / (base::rowSums(cov_mat[, control_idx]) + 1e-100)

    #####################################

    result = do.call(rbind, parallel::mclapply(seq_along(gr), function(locus_idx){

        ### Deal with local information (or not)
        if(local_window_size != 0) {
            # NOTE: It is much faster to work with subsets of the result of start()
            # than it is to work with subsets of GRanges.

            # Get the indices which are within the local_window_size, but also limit to 5 CpGs on either side
            # NOTE, local information is only used with cytosine/CpG resolution data so start() is valid.
            # If regions were allowed, we would have to pay attention to which side we're on and use start()/end()
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

        #####################################

        ### Compute the degrees of freedom for the locus
        if(all(disp_groups)) {
            df_subtract = 2
        } else {
            df_subtract = 1
        }
        df = pmax(rowSums(local_cov[, disp_groups_idx, drop = FALSE] > 0) - df_subtract, 0)
        # Compute the degrees of freedom to be used in the test for differential methylation
        df = sum(df * local_weights)

        #####################################

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

            #####################################

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

            #####################################

            ### log Likelihood ratio calculation
            log_lik_ratio =
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

            #####################################

            locus_data = c(
                disp_est = disp_est,
                log_lik_ratio = log_lik_ratio,
                meth_control = group_meth_est[1]*100,
                meth_case = group_meth_est[2]*100,
                meth_all = group_meth_est[3]*100,
                df = df + 2)
        } else {
            # Not enough degrees of freedom, return NAs, these will be removed
            # with a message to the user with how many
            locus_data = c(
                disp_est = NA,
                log_lik_ratio = NA,
                meth_control = NA,
                meth_case = NA,
                meth_all = NA,
                df = df)
        }

        return(locus_data)
    }, mc.cores = n_cores))

    #####################################

    # Build GRanges version of result
    result_gr = gr
    mcols(result_gr) = result

    # Calculate pvalue
    if(t_approx) {
         result_gr$pvalue = stats::pt(-sqrt(pmax(result_gr$log_lik_ratio, 0)), result_gr$df) * 2
    } else {
         result_gr$pvalue = stats::pchisq(pmax(result_gr$log_lik_ratio, 0), 1, lower.tail = FALSE)
    }

    # Calculate meth_diff and set very small differences to 0
    result_gr$meth_diff = (result_gr$meth_case - result_gr$meth_control)
    result_gr$meth_diff[abs(result_gr$meth_diff) < 0.01] = 0

    # Correct for multiple testing
    result_gr$fdr = stats::p.adjust(result_gr$pvalue, method = 'BH')

    # Assign direction of hypermethylation (NOTE, this is not "significance")
    result_gr$direction = ifelse(result_gr$meth_diff >= 0, case, control)

    #####################################

    # Order output columns and attach to GRanges
    col_order = c(
        'meth_case',
        'meth_control',
        'meth_diff',
        'direction',
        'pvalue',
        'fdr',
        'disp_est',
        'log_lik_ratio',
        'df'
    )
    mcols(result_gr) = mcols(result_gr)[, col_order]

    #####################################

    # Check for NA results and indicate how many loci were dropped because of
    # a lack of available degrees of freedom
    insufficient_df = result_gr$df == 1

    if(sum(insufficient_df) > 0) {
        result_gr = result_gr[!insufficient_df]
        message(sprintf('%s loci were dropped due to insufficient degrees of freedom (df = 1).', sum(insufficient_df)))
    }

    return(result_gr)
}
