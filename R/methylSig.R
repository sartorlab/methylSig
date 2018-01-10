#' Differential methylation analysis using binomial model
#'
#' This function calculates differential methylation statistics using a binomial-based approach. See `Warning' message below.
#'
#' This function uses a binomial-based model to calculate differential methylation statistics. It is nearly identical to the \code{methylKit::calculateDiffMeth} function in the \code{methylKit} R package except that only the likelihood ratio test and \code{p.adjust()} with \code{method=``BH''} are used to calculate significance levels. It is significantly faster than \code{methylKit::calculateDiffMeth} function.
#'
#' @param meth A \code{methylSigData-class} object to calculate differential methylation statistics. It can be obtained using \code{\link{methylSigReadData}}.
#' @param groups A vector of two numbers specify two groups to compare. See \code{treatment} argument of \code{\link{methylSigReadData}} function. Default is \code{c(Treatment=1,Control=0)}.
#' @param min.per.group A vector with two numbers that specify the minimum numbers of samples required to be qualify as defferentially methylation region. If it is a single number, both groups will use it as the minimum requried number of samples. Default is \code{c(3,3)}.
#'
#' @return A \code{methylSigDiff-class} object that contains the differential methylation statistics and chromosomal locations. \code{p.adjust} with \code{method="BH"} option is used for P-value correction.
#'
#' @section Warning: This function does not take into account the variability among samples in each group being compared.
#'
#' @seealso \code{\link{methylSigCalc}}
#'
#' @examples
#' data(sampleData)
#' myDiff = binomialDiffCalc(meth)
#'
#' @keywords differentialMethylation
#'
#' @export
binomialDiffCalc <- function(meth, groups=c("Treatment"=1,"Control"=0), min.per.group=c(3,3)) {
    treatment = slot(meth, "treatment")

    group1 = which(treatment == groups[1])
    group2 = which(treatment == groups[2])

    if(length(min.per.group) == 1) {
        min.per.group = c(min.per.group,min.per.group)
    }

    vlist = rowSums(meth@data.coverage[,group1]>0,na.rm = TRUE) >= min.per.group[1] &
                   rowSums(meth@data.coverage[,group2]>0, na.rm=TRUE) >= min.per.group[2]
    treads1 = rowSums(meth@data.numTs[vlist,group1], na.rm = TRUE)
    treads2 = rowSums(meth@data.numTs[vlist,group2], na.rm = TRUE)
    creads1 = rowSums(meth@data.numCs[vlist,group1], na.rm = TRUE)
    creads2 = rowSums(meth@data.numCs[vlist,group2], na.rm = TRUE)

    logLikRatio <- 2*(creads1*log(creads1/(treads1+creads1)+1e-100)
                      + treads1*log(treads1/(treads1+creads1)+1e-100)
                      + creads2*log(creads2/(treads2+creads2)+1e-100)
                      + treads2*log(treads2/(treads2+creads2)+1e-100)
                      - (creads1 + creads2)*log((creads1+creads2)/(creads1 + creads2 + treads1 + treads2)+1e-100)
                      - (treads1 + treads2)*log((treads1+treads2)/(creads1 + creads2 + treads1 + treads2)+1e-100)
                     )

    pvalue = pchisq(logLikRatio, 1, lower.tail=FALSE)

    results=cbind(pvalue,p.adjust(pvalue,method ="BH"),(creads1/(creads1+treads1)-creads2/(creads2+treads2))*100,logLikRatio,
                  creads1/(creads1+treads1)*100, creads2/(creads2+treads2)*100)
    colnames(results) = c("pvalue","qvalue", "meth.diff","logLikRatio", paste("mu", groups,sep=""))
    methylSig.newDiff(meth@data.ids[vlist], meth@data.chr[vlist], meth@data.start[vlist],meth@data.end[vlist],
                              meth@data.strand[vlist], results, sample.ids=meth@sample.ids[c(group1,group2)],
                              sample.filenames=meth@sample.filenames[c(group1,group2)],
                              treatment=meth@treatment[c(group1,group2)], destranded=meth@destranded,
                              resolution=meth@resolution,
                              options=paste("min.per.group=c(",min.per.group[1],",",min.per.group[2],") & Total: ", NROW(results), sep=""),
                              data.options = meth@options)
}

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
#' @param meth A \code{methylSigData-class} object to calculate differential methylation statistics. It can be obtained using `methylSigReadData'.
#' @param comparison The name of the column in \code{pData(meth)} to use for the comparisons.
#' @param dispersion A value indicating which group or groups are used to estimate the variance (dispersion). If groups are defined as c(Treatment=1,Control=0), dispersion can take values "Treatment", "Control", 1, 0 or "both". Default is "both".
#' @param local.info A as.logical value indicating whether to use local information to improve dispersion parameter estimation. Default is FALSE.
#' @param local.winsize A number to specify the window size in basepairs for local dispersion estimation. The dispersion (variance) parameter of the groups at the particular location LOC is calculated using information from LOC - local.winsize to LOC + local.winsize. This argument is only activated when local.info=TRUE. Default is 200.
#' @param min.per.group A vector with two numbers that specify the minimum numbers of samples required to be qualify as defferentially methylated region.  If it is a single number, both groups will use it as the minimum requried number of samples. Default is c(3,3).
#' @param weightFunc A weight kernel function. The input of this function is from -1 to 1. The default is the tri-weight kernel function defined as function(u) = (1-u^2)^3. Function value and range of parameter for weight function should be from 0 to 1.
#' @param T.approx A as.logical value indicating whether to use squared t approximation for the likelihood ratio statistics. Chi-square approximation (T.approx = FALSE) is recommended when the sample size is large.  Default is TRUE.
#' @param num.cores An integer denoting how many cores should be used for differential methylation calculations (only can be used in machines with multiple cores).
#'
#' @return \code{methylSigDiff-class} object containing the differential methylation statistics and locations. p.adjust with method="BH" option is used for P-value correction.
#'
#' @seealso \code{\link{methylSigPlot}}, \code{\link{methylSigReadData}}
#'
#' @examples
#' data(sampleData)
#'
#' myDiffSig = methylSigCalc (meth)
#'
#' ### calculate differential methylation statistics using
#' ### treatment group 0 to evaluate dispersion parameter.
#' ### Also use local information to improve variance estimation.
#'
#' myDiffSig = methylSigCalc (meth, dispersion=0, local.info=TRUE,
#'                            min.per.group=4)
#'
#' @keywords differentialMethylation
#'
#' @export
methylSigCalc = function(meth, comparison = NA, dispersion="both",
         local.info=FALSE, local.winsize=200,
         min.per.group=c(3,3), weightFunc=methylSig_weightFunc, T.approx = TRUE,
         num.cores = 1) {

    if(!local.info) {
        local.winsize = 0
    }

    min.disp = 1e-6
    min.InvDisp = 0.001
    max.InvDisp = max(1/max(min.disp, 1e-6), min.InvDisp)

    minMu = 0
    maxMu = 1

    # Get the group labels, THIS ASSUMES CORRECT REFERENCE LEVEL SET
    pdata = pData(meth)
    group2 = levels(pdata[, comparison])[2]
    group1 = levels(pdata[, comparison])[1]

    # Determine which rows of pData belong to which group
    # / which columns of Cov and M matrices belong to which group
    group2_idx = which(pdata[,comparison] == group2)
    group1_idx = which(pdata[,comparison] == group1)

    # Determine which sample column indexes to use for dispersion calculation
    if(dispersion == 'both') {
        local_idx = c(group2_idx, group1_idx)
    } else if (dispersion == group2) {
        local_idx = group2_idx
    } else if (dispersion == group1) {
        local_idx = group1_idx
    } else {
        stop('"dispersion" should be one of "both", the name of group2, or the name of group1')
    }

    meth_gr = GenomicRanges::granges(meth)

    all_cov = as.matrix(bsseq::getCoverage(meth, type = 'Cov'))
    all_meth = as.matrix(bsseq::getCoverage(meth, type = 'M'))

    muEst = matrix(0, ncol = ncol(meth), nrow = nrow(meth))
    muEst[, group2_idx] = base::rowSums(all_meth[, group2_idx]) / (base::rowSums(all_cov[, group2_idx]) + 1e-100)
    muEst[, group1_idx] = base::rowSums(all_meth[, group1_idx]) / (base::rowSums(all_cov[, group1_idx]) + 1e-100)

    # Determine which loci satisfy min.per.group
    valid_idx = which(
        base::rowSums(all_cov[, group2_idx] > 0) >= min.per.group[2] & base::rowSums(all_cov[, group1_idx] > 0) >= min.per.group[1]
    )

    # Go through each index for a valid locus
    results = do.call(rbind, mclapply(valid_idx, function(idx){

        if(local.winsize != 0) {
            # Use local information
            locus = meth_gr[idx]

            # Setup GRanges where we expand from each CpG by the winsize
            # This is a GRanges object with only 1 range
            local_win_gr = GenomicRanges::flank(locus, width = local.winsize, start = TRUE, both = TRUE)

            # Look for overlaps in the original meth for each of the above two
            # NOTE: You can use meth (BSseq) and it will access the GRanges for GenomicRanges::findOverlaps()
            local_overlaps = GenomicRanges::findOverlaps(query = meth_gr, subject = local_win_gr, ignore.strand = TRUE)

            query_idx = S4Vectors::queryHits(local_overlaps)

            if(length(query_idx) == 1) {
                # Do not use local information
                query_idx = idx
                local_weights = 1

                # Collect Cov and M matrices for all the loci in the window
                # These are delayed matrices. Rows are loci and columns are samples
                local_cov = matrix(all_cov[query_idx, ], nrow = 1)
                local_meth = matrix(all_meth[query_idx, ], nrow = 1)

                # Collect the correct rows of muEst
                local_muEst = matrix(muEst[query_idx, ], nrow = 1)
            } else {
                # Normalize the locations of the loci in the window to be in the domain of
                # the weight function [-1, 1]
                # locus is GRanges with a single range
                # local_loci and meth_loci are GRanges with number of ranges equal to overlaps
                # based on the window size
                local_loci = meth_gr[query_idx]

                # Each is a vector of input values to the weight function
                # We need to scale the loci in the window onto the interval [-1, 1] because
                # that is the domain of the weightFuncction.
                local_loci_norm = (start(local_loci) - start(locus)) / (local.winsize + 1)

                # Calculate the weights
                # Each is a vector of values of the weight function.
                local_weights = weightFunc(local_loci_norm)

                # Collect Cov and M matrices for all the loci in the window
                # These are delayed matrices. Rows are loci and columns are samples
                local_cov = all_cov[query_idx, ]
                local_meth = all_meth[query_idx, ]

                # Collect the correct rows of muEst
                local_muEst = muEst[query_idx, ]
            }
        } else {
            # Do not use local information
            query_idx = idx
            local_weights = 1

            # Collect Cov and M matrices for all the loci in the window
            # These are delayed matrices. Rows are loci and columns are samples
            local_cov = matrix(all_cov[query_idx, ], nrow = 1)
            local_meth = matrix(all_meth[query_idx, ], nrow = 1)

            # Collect the correct rows of muEst
            local_muEst = matrix(muEst[query_idx, ], nrow = 1)
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
        df = pmax(rowSums(local_cov[, local_idx, drop = FALSE] > 0) - df_subtract, 0)
        # Compute the degrees of freedom to be used in the test for differential methylation
        df = sum(df * local_weights)

        if(df > 1) {
            ### Common dispersion calculation
            # This returns a singleton numeric
            if(methylSig_derivativePhi(
                phi = max.InvDisp,
                lCreads = local_creads[, local_idx, drop = FALSE],
                lTreads = local_treads[, local_idx, drop = FALSE],
                mu = local_muEst[, local_idx, drop = FALSE],
                weight = local_weights) >= 0) {

                # Describe
                phiCommonEst = max.InvDisp
            } else if(methylSig_derivativePhi(
                phi = min.InvDisp,
                lCreads = local_creads[, local_idx, drop = FALSE],
                lTreads = local_treads[, local_idx, drop = FALSE],
                mu = local_muEst[, local_idx, drop = FALSE],
                weight = local_weights) <= 0){

                # Describe
                phiCommonEst = min.InvDisp
            } else {
                # Describe
                phiCommonEst = stats::uniroot(
                    f = methylSig_derivativePhi,
                    interval = c(min.InvDisp, max.InvDisp),
                    local_creads[, local_idx, drop = FALSE],
                    local_treads[, local_idx, drop = FALSE],
                    local_muEst[, local_idx, drop = FALSE],
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

            locus_data = data.frame(
                phiCommonEst = phiCommonEst,
                logLikRatio = logLikRatio,
                muEstC_group1 = muEstC[1]*100,
                muEstC_group2 = muEstC[2]*100,
                muEstC_group12 = muEstC[3]*100,
                df = df + 2,
                stringsAsFactors = F)
        } else {
            # Not enough degrees of freedom, return NAs
            locus_data = data.frame(
                phiCommonEst = NA,
                logLikRatio = NA,
                muEstC_group1 = NA,
                muEstC_group2 = NA,
                muEstC_group12 = NA,
                df = NA,
                stringsAsFactors = F)
        }

        return(locus_data)
    }, mc.cores = num.cores))

    results_gr = GenomicRanges::granges(meth)[valid_idx]
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

    return(results_gr)
}
