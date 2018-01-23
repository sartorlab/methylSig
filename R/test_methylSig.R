# Called by methylSigCalc
methylSig_weightFunc <- function(u) (1-u^2)^3

# Called by methylSigCalc
methylSig_derivativePhi <- function(phi, lCreads, lTreads, mu, weight) {
    derivative = 0
    indicator_c = lCreads > 0
    indicator_t = lTreads > 0
    indicator_ct = lCreads + lTreads > 0

    if(nrow(lCreads) == 1) {
        derivative =
            sum( indicator_c * ( mu * (digamma((mu * phi) + lCreads + 1e-100) - digamma(mu * phi + 1e-100)) ) ) +
            sum( indicator_t * ((1 - mu) * (digamma( ((1 - mu) * phi) + lTreads + 1e-100) - digamma( ((1-mu) * phi) + 1e-100))) ) -
            sum( indicator_ct * (digamma(phi + lCreads + lTreads + 1e-100) - digamma(phi)) )
    } else {
        for(g in 1:ncol(lCreads)) {
            derivative = derivative +
                sum( indicator_c[,g] * (weight * mu[,g] * (digamma(mu[,g] * phi + lCreads[,g] + 1e-100) - digamma(mu[,g] * phi + 1e-100))) ) +
                sum( indicator_t[,g] * (weight * (1 - mu[,g]) * (digamma((1 - mu[,g]) * phi + lTreads[,g] + 1e-100) - digamma((1 - mu[,g]) * phi + 1e-100))) ) -
                sum( indicator_ct[,g] * (weight * (digamma(phi + lCreads[,g] + lTreads[,g] + 1e-100) - digamma(phi))) )
        }
    }

    derivative
}

# Called by methylSigCalc
methylSig_derivativeMu <- function(mu, lCreads, lTreads, phi, weight) {
    derivative = 0
    indicator_c = lCreads > 0
    indicator_t = lTreads > 0

    if(nrow(lCreads) == 1) {
        derivative =
            sum( indicator_c * (digamma(mu * phi + lCreads + 1e-100) - digamma(mu * phi + 1e-100)) ) -
            sum( indicator_t * (digamma((1 - mu) * phi + lTreads + 1e-100) - digamma((1 - mu) * phi + 1e-100)) )
    } else {
        for(g in 1:ncol(lCreads)) {
            derivative = derivative +
                sum( indicator_c[,g] * (weight * (digamma(mu * phi + lCreads[,g]+ 1e-100) - digamma(mu * phi + 1e-100))) ) -
                sum( indicator_t[,g] * (weight * (digamma((1 - mu) * phi + lTreads[,g] + 1e-100) - digamma((1 - mu) * phi + 1e-100))) )
        }
    }

    derivative
}

# Called by methylSig_dataProcess
methylSig_logLik  <- function(mu, phi, lCreads, lTreads, weight) {
    llik = 0
    indicator_c = lCreads > 0
    indicator_t = lTreads > 0

    if(nrow(lCreads) == 1) {
        llik = llik +
            sum( indicator_c * (lgamma(mu * phi + lCreads + 1e-100) - lgamma(mu * phi + 1e-100)) ) +
            sum( indicator_t * (lgamma((1 - mu) * phi + lTreads + 1e-100) - lgamma((1 - mu) * phi + 1e-100)) )
    } else {
        for(g in 1:ncol(lCreads)) {
            llik = llik +
                sum( indicator_c[,g] * (weight * (lgamma(mu * phi + lCreads[,g] + 1e-100) - lgamma(mu * phi + 1e-100))) ) +
                sum( indicator_t[,g] * (weight * (lgamma((1 - mu) * phi + lTreads[,g] + 1e-100) - lgamma((1 - mu) + 1e-100))) )
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
#' @param meth A \code{BSseq-class} object to calculate differential methylation statistics. See \code{methylSigReadData} for how to read in methylation data.
#' @param comparison The name of the column in \code{pData(meth)} to use for the comparisons.
#' @param dispersion One of \code{both}, or either group name. Indicates which set of samples to use to estimate the dispersion parameter. Default is \code{both}.
#' @param local.info A \code{logical} value indicating whether to use local information to improve mean and dispersion parameter estimations. Default is \code{FALSE}.
#' @param local.winsize An \code{integer} to specify the distance upstream and downstream of a location to include local information for the mean and dispersion parameter estimations. NOTE: An additional constraint is placed whereby a maximum of 5 loci upstream and downstream of the locus of interest are used. Default is \code{200}.
#' @param min.per.group A vector with two numbers specifying the minimum number of samples required to perform the test for differential methylation. If it is a single number, both groups will use it as the minimum requried number of samples. Default is \code{c(3,3)}.
#' @param weightFunc A weight kernel function. The input of this function is from -1 to 1. The default is the tri-weight kernel function defined as \code{function(u) = (1-u^2)^3}. Function value and range of parameter for weight function should be from 0 to 1.
#' @param T.approx A \code{logical} value indicating whether to use squared t approximation for the likelihood ratio statistics. Chi-square approximation (\code{T.approx = FALSE}) is recommended when the sample size is large.  Default is \code{TRUE}.
#' @param num.cores An integer denoting how many cores should be used for differential methylation calculations.
#'
#' @return A \code{GRanges} object containing the following \code{mcols}:
#' \describe{
#'   \item{phiCommonEst}{ The dispersion estimate. }
#'   \item{logLikRatio}{ The log likelihood ratio. }
#'   \item{df}{ Degrees of freedom used when \code{T.approx = TRUE}. }
#'   \item{muEstC_group1}{ Methylation estimate for group1. Groups correspond to the levels in the column used for the comparison in \code{pdata(meth)}. }
#'   \item{muEstC_group2}{ Methylation estimate for group2. }
#'   \item{meth.diff}{ The difference \code{muEstC_group2 - muEstC_group1}. }
#'   \item{hyper_direction}{ The group for which the CpG/region is hyper-methylated. Groups correspond to the levels in the column used for the comparison in \code{pdata(meth)}. }
#'   \item{pvalue}{ The p-value from the t-test (\code{T.approx = TRUE}) or the Chi-Square test (\code{T.approx = FALSE}). }
#'   \item{fdr}{ The Benjamini-Hochberg adjusted p-values using \code{p.adjust(method = 'BH')}. }
#' }
#'
#' @seealso \code{\link{methylSigPlot}}, \code{\link{methylSigReadData}}
#'
#' @examples
#' data(data, package = 'methylSig')
#'
#' result = methylSigCalc(
#'     meth = data,
#'     comparison = 'DR_vs_DS',
#'     dispersion = 'both',
#'     local.info = FALSE,
#'     local.winsize = 200,
#'     min.per.group = c(3,3),
#'     weightFunc = methylSig_weightFunc,
#'     T.approx = TRUE,
#'     num.cores = 1)
#'
#' @keywords differentialMethylation
#'
#' @export
methylSigCalc = function(
    meth,
    comparison = NA,
    dispersion="both",
    local.info=FALSE, local.winsize=200,
    min.per.group=c(3,3),
    weightFunc=methylSig_weightFunc,
    T.approx = TRUE,
    num.cores = 1) {

    #####################################
    # Constants
    if(!local.info) {
        local.winsize = 0
    }

    if(length(min.per.group) == 1) {
        min.per.group = c(min.per.group,min.per.group)
    }

    min.disp = 1e-6
    min.InvDisp = 0.001
    max.InvDisp = max(1/max(min.disp, 1e-6), min.InvDisp)

    minMu = 0
    maxMu = 1

    #####################################
    # Get the group labels, NOTE: THIS ASSUMES CORRECT REFERENCE LEVEL SET
    pdata = pData(meth)
    group2 = levels(pdata[, comparison])[2]
    group1 = levels(pdata[, comparison])[1]

    # Determine which rows of pData belong to which group
    # / which columns of Cov and M matrices belong to which group
    group2_idx = which(pdata[,comparison] == group2)
    group1_idx = which(pdata[,comparison] == group1)

    # Determine which sample column indexes to use for dispersion calculation
    if(dispersion == 'both') {
        disp_groups_idx = c(group2_idx, group1_idx)
    } else if (dispersion == group2) {
        disp_groups_idx = group2_idx
    } else if (dispersion == group1) {
        disp_groups_idx = group1_idx
    } else {
        stop('"dispersion" should be one of "both", the name of group2, or the name of group1')
    }

    #####################################
    # Determine which sites are valid to test according to min.per.group
    all_cov = as.matrix(bsseq::getCoverage(meth, type = 'Cov'))
    all_meth = as.matrix(bsseq::getCoverage(meth, type = 'M'))

    # Estimate mu per locus within each group. The same value is used for all samples within the same group.
    muEst = matrix(0, ncol = ncol(meth), nrow = nrow(meth))
    muEst[, group2_idx] = base::rowSums(all_meth[, group2_idx]) / (base::rowSums(all_cov[, group2_idx]) + 1e-100)
    muEst[, group1_idx] = base::rowSums(all_meth[, group1_idx]) / (base::rowSums(all_cov[, group1_idx]) + 1e-100)

    # Determine which loci satisfy min.per.group
    valid_idx = which(
        base::rowSums(all_cov[, group2_idx] > 0) >= min.per.group[2] & base::rowSums(all_cov[, group1_idx] > 0) >= min.per.group[1]
    )

    #####################################
    # Resize all_cov, all_meth, and muEst to valid_idx
    # Extract the granges of the meth BSseq object
    # These are all used within the mclapply below
    all_cov = all_cov[valid_idx,]
    all_meth = all_meth[valid_idx,]
    muEst = muEst[valid_idx,]
    meth_gr = granges(meth)[valid_idx,]

    num_loci = length(valid_idx)

    #####################################
    # Go through each valid locus
    results = do.call(rbind, mclapply(seq_along(valid_idx), function(locus_idx){

        if(local.winsize != 0) {
            # NOTE: It is much faster to work with subsets of the result of start()
            # than it is to work with subsets of GRanges.

            # Get the indices which are within the local.winsize, but also limit to 5 CpGs on either side
            local_loci_idx = intersect(
                which(abs(start(meth_gr)[locus_idx] - start(meth_gr)) < local.winsize),
                max(1, locus_idx - 5):min(num_loci, locus_idx + 5))

            if(length(local_loci_idx) == 1) {
                # Do not use local information
                local_loci_idx = locus_idx
                local_weights = 1

                # Collect Cov and M matrices for all the loci in the window
                # Rows are loci and columns are samples
                local_cov = matrix(all_cov[local_loci_idx, ], nrow = 1)
                local_meth = matrix(all_meth[local_loci_idx, ], nrow = 1)

                # Collect the correct rows of muEst
                local_muEst = matrix(muEst[local_loci_idx, ], nrow = 1)
            } else {
                # We need to scale the loci in the window onto the interval [-1, 1] because
                # that is the domain of the weightFunc.
                # This is a vector of the distances of the local loci to the loci of interest (domain)
                local_loci_norm = (start(meth_gr)[local_loci_idx] - start(meth_gr)[locus_idx]) / (local.winsize + 1)

                # Calculate the weights
                # Each is a vector of values of the weight function (range)
                local_weights = weightFunc(local_loci_norm)

                # Collect Cov and M matrices for all the loci in the window
                # Rows are loci and columns are samples
                local_cov = all_cov[local_loci_idx, ]
                local_meth = all_meth[local_loci_idx, ]

                # Collect the correct rows of muEst
                local_muEst = muEst[local_loci_idx, ]
            }
        } else {
            # Do not use local information
            local_loci_idx = locus_idx
            local_weights = 1

            # Collect Cov and M matrices for all the loci in the window
            # Rows are loci and columns are samples
            local_cov = matrix(all_cov[local_loci_idx, ], nrow = 1)
            local_meth = matrix(all_meth[local_loci_idx, ], nrow = 1)

            # Collect the correct rows of muEst
            local_muEst = matrix(muEst[local_loci_idx, ], nrow = 1)
        }

        # Convert to old methylSig notion of C reads (methylated) and T reads (unmethylated)
        # Then we can reuse the derivative and log likelihood functions Yongseok implemented.
        # These are matrices. Rows are loci and columns are samples
        local_creads = local_meth
        local_treads = local_cov - local_meth

        # Compute the degrees of freedom for the locus
        if(dispersion == 'both') {
            df_subtract = 2
        } else {
            df_subtract = 1
        }
        df = pmax(rowSums(local_cov[, disp_groups_idx, drop = FALSE] > 0) - df_subtract, 0)
        # Compute the degrees of freedom to be used in the test for differential methylation
        df = sum(df * local_weights)

        if(df > 1) {
            ### Common dispersion calculation
            # This returns a singleton numeric
            if(methylSig_derivativePhi(
                phi = max.InvDisp,
                lCreads = local_creads[, disp_groups_idx, drop = FALSE],
                lTreads = local_treads[, disp_groups_idx, drop = FALSE],
                mu = local_muEst[, disp_groups_idx, drop = FALSE],
                weight = local_weights) >= 0) {

                phiCommonEst = max.InvDisp
            } else if(methylSig_derivativePhi(
                phi = min.InvDisp,
                lCreads = local_creads[, disp_groups_idx, drop = FALSE],
                lTreads = local_treads[, disp_groups_idx, drop = FALSE],
                mu = local_muEst[, disp_groups_idx, drop = FALSE],
                weight = local_weights) <= 0){

                phiCommonEst = min.InvDisp
            } else {
                phiCommonEst = stats::uniroot(
                    f = methylSig_derivativePhi,
                    interval = c(min.InvDisp, max.InvDisp),
                    local_creads[, disp_groups_idx, drop = FALSE],
                    local_treads[, disp_groups_idx, drop = FALSE],
                    local_muEst[, disp_groups_idx, drop = FALSE],
                    local_weights)$root
            }

            ### Common group means calculation
            # This returns a numeric vector (group1, group2, group1 + group2) with the mu estimate
            muEstC_groups = list(group1_idx, group2_idx, c(group1_idx, group2_idx))
            muEstC = rep(0, length(muEstC_groups))
            for(group_idx in seq_along(muEstC_groups)) {
                if(sum(local_creads[, muEstC_groups[[group_idx]], drop = FALSE]) == 0) {
                    # If there are no local C reads, methylation is 0
                    muEstC[group_idx] = 0
                } else if (sum(local_treads[, muEstC_groups[[group_idx]], drop = FALSE]) == 0) {
                    # If there are no local T reads, methylation is 1
                    muEstC[group_idx] = 1
                } else {
                    # Otherwise, do something fancier
                    muEstC[group_idx] = stats::uniroot(
                        f = methylSig_derivativeMu,
                        interval = c(minMu, maxMu),
                        local_creads[, muEstC_groups[[group_idx]], drop = FALSE],
                        local_treads[, muEstC_groups[[group_idx]], drop = FALSE],
                        phiCommonEst,
                        local_weights)$root
                }
            }

            ### log Likelihood ratio calculation
            logLikRatio =
                methylSig_logLik(
                    mu = muEstC[1],
                    phi = phiCommonEst,
                    lCreads = local_creads[, group1_idx, drop = FALSE],
                    lTreads = local_treads[, group1_idx, drop = FALSE],
                    weight = local_weights) +
                methylSig_logLik(
                    mu = muEstC[2],
                    phi = phiCommonEst,
                    lCreads = local_creads[, group2_idx, drop = FALSE],
                    lTreads = local_treads[, group2_idx, drop = FALSE],
                    weight = local_weights) -
                methylSig_logLik(
                    mu = muEstC[3],
                    phi = phiCommonEst,
                    lCreads = local_creads[, c(group1_idx, group2_idx), drop = FALSE],
                    lTreads = local_treads[, c(group1_idx, group2_idx), drop = FALSE],
                    weight = local_weights)

            locus_data = c(
                phiCommonEst = phiCommonEst,
                logLikRatio = logLikRatio,
                muEstC_group1 = muEstC[1]*100,
                muEstC_group2 = muEstC[2]*100,
                muEstC_group12 = muEstC[3]*100,
                df = df + 2)
        } else {
            # Not enough degrees of freedom, return NAs
            locus_data = c(
                phiCommonEst = NA,
                logLikRatio = NA,
                muEstC_group1 = NA,
                muEstC_group2 = NA,
                muEstC_group12 = NA,
                df = NA)
        }

        return(locus_data)
    }, mc.cores = num.cores))

    results_gr = meth_gr
    mcols(results_gr) = results

    if(T.approx) {
         results_gr$pvalue = stats::pt(-sqrt(pmax(results_gr$logLikRatio, 0)), results_gr$df) * 2
    } else {
         results_gr$pvalue = stats::pchisq(pmax(results_gr$logLikRatio, 0), 1, lower.tail = F)
    }

    # Set any methylation difference less than 0.01 to 0
    results_gr$meth.diff = (results_gr$muEstC_group2 - results_gr$muEstC_group1)
    results_gr$meth.diff[abs(results_gr$meth.diff) < 0.01] = 0
    results_gr$meth.diff = as.numeric(results_gr$meth.diff)

    results_gr$fdr = p.adjust(results_gr$pvalue, method = 'BH')

    results_gr$hyper_direction = ifelse(results_gr$meth.diff > 0, group2, group1)

    mcols(results_gr) = mcols(results_gr)[, c('phiCommonEst', 'logLikRatio', 'df', 'muEstC_group1', 'muEstC_group2', 'meth.diff', 'hyper_direction', 'pvalue', 'fdr')]

    return(results_gr)
}
