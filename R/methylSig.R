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

methylSig_weightFunc <- function(u) (1-u^2)^3

# Called by methylSig_dataProcess
methylSig_derivativePhi <- function(phi, lCreads, lTreads, mu, weight) {
    derivative = 0
    if(nrow(lCreads) == 1) {
        ### Only one location, weight does not matter
            vlist <- which(as.logical(lCreads > 0))
            if(length(vlist) > 0)
                derivative = derivative + sum( mu[vlist] * (digamma((mu[vlist] * phi) + lCreads[vlist] + 1e-100) - digamma(mu[vlist] * phi + 1e-100)) )

            vlist <- which(as.logical(lTreads > 0))
            if(length(vlist) > 0)
                derivative = derivative + sum( (1 - mu[vlist]) * (digamma( ((1 - mu[vlist]) * phi) + lTreads[vlist] + 1e-100) - digamma( ((1-mu[vlist]) * phi) + 1e-100)))

            vlist <- which(as.logical(lCreads + lTreads) > 0)
            if(length(vlist) > 0)
                derivative = derivative - sum( digamma(phi + lCreads[vlist] + lTreads[vlist] + 1e-100) - digamma(phi))
    } else {
        for(g in 1:ncol(lCreads)) {
            vlist <- which(as.logical(lCreads[, g] > 0))
            if(length(vlist) > 0)
                derivative = derivative + sum( weight[vlist] * mu[vlist, g] * (digamma(mu[vlist, g] * phi + lCreads[vlist, g] + 1e-100) - digamma(mu[vlist, g] * phi + 1e-100)) )

            vlist <- which(as.logical(lTreads[vlist, g] > 0))
            if(length(vlist) > 0)
                derivative = derivative + sum( weight[vlist] * (1 - mu[vlist, g]) * (digamma((1 - mu[vlist, g]) * phi + lTreads[vlist, g] + 1e-100) - digamma((1 - mu[vlist, g]) * phi + 1e-100)) )

            vlist <- which(as.logical((lCreads[, g] + lTreads[, g]) > 0))
            if(length(vlist) > 0)
                derivative = derivative - sum( weight[vlist] * (digamma(phi + lCreads[vlist, g] + lTreads[vlist, g] + 1e-100) - digamma(phi)) )
        }
    }

    derivative
}

# Called by methylSig_dataProcess
methylSig_derivativeMu <- function(mu, lCreads, lTreads, phi, weight) {
    derivative <- 0
    if(nrow(lCreads) == 1) {
        vlist <- which(as.logical(lCreads > 0))
        if(length(vlist) > 0)
            derivative = derivative + sum(digamma(mu * phi + lCreads[vlist] + 1e-100) - digamma(mu * phi + 1e-100))
        vlist <- which(as.logical(lTreads > 0))
        if(length(vlist) > 0)
            derivative = derivative - sum(digamma((1 - mu) * phi + lTreads[vlist] + 1e-100) - digamma((1 - mu) * phi + 1e-100))
    } else {
        for(g in 1:ncol(lCreads)) {
            vlist <- which(as.logical(lCreads[,g] > 0))
            if(length(vlist) > 0)
                derivative = derivative + sum(weight[vlist] * (digamma(mu * phi + lCreads[vlist, g]+ 1e-100) - digamma(mu * phi + 1e-100)))
            vlist <- which(as.logical(lTreads[, g] > 0))
            if(length(vlist) > 0)
                derivative = derivative - sum(weight[vlist] * (digamma((1 - mu) * phi + lTreads[vlist, g] + 1e-100) - digamma((1 - mu) * phi + 1e-100)))
        }
    }

    derivative
}

# Called by methylSig_dataProcess
methylSig_logLik  <- function(mu, phi, lCreads, lTreads, weight) {
    llik = 0
    if(nrow(lCreads) == 1) {
        ###### single location, weight should be 1
        for(g in 1:ncol(lCreads)) {
            if(lCreads[g] > 0) {
                llik = llik + as.numeric(lgamma(mu * phi + lCreads[g]) - lgamma(mu * phi + 1e-100))
            }
            if(lTreads[g] > 0) {
                llik = llik + as.numeric(lgamma((1 - mu) * phi + lTreads[g]) - lgamma((1 - mu) * phi + 1e-100))
            }
        }
    } else {
        for(g in 1:ncol(lCreads)) {
            vlist <- which(as.logical(lCreads[,g] > 0))
            if(length(vlist) > 0) {
                for(i in vlist) {
                    llik = llik + as.numeric(weight[i] * (lgamma(mu * phi + lCreads[i, g]) - lgamma(mu * phi + 1e-100)))
                }
            }
            vlist <- which(as.logical(lTreads[, g] > 0))
            if(length(vlist) > 0) {
                for(i in vlist) {
                    llik = llik + as.numeric(weight[i] * (lgamma((1 - mu) * phi + lTreads[i, g]) - lgamma((1 - mu) + 1e-100)))
                }
            }
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
#' @param local.disp A as.logical value indicating whether to use local information to improve dispersion parameter estimation. Default is FALSE.
#' @param winsize.disp A number to specify the window size in basepairs for local dispersion estimation. The dispersion (variance) parameter of the groups at the particular location LOC is calculated using information from LOC - winsize.disp to LOC + winsize.disp. This argument is only activated when local.disp=TRUE. Default is 200.
#' @param local.meth A as.logical value indicating whether to use local information to improve methylation level estimation. Default is FALSE.
#' @param winsize.meth a number to specify the window size in basepairs for local methylation level estimation. The group methylation level at the particular location LOC is calculated using information from LOC - winsize.meth to LOC + winsize.meth. This argument is only activated when local.meth=TRUE. Default is 200.
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
#' myDiffSig = methylSigCalc (meth, dispersion=0, local.disp=TRUE,
#'                            min.per.group=4)
#'
#' @keywords differentialMethylation
#'
#' @export
methylSigCalc = function(meth, comparison = NA, dispersion="both",
         local.disp=FALSE, winsize.disp=200,
         local.meth=FALSE, winsize.meth=200,
         min.per.group=c(3,3), weightFunc=methylSig_weightFunc, T.approx = TRUE,
         num.cores = 1) {

    if(!local.disp) {
        winsize.disp = 0
    }
    if(!local.meth) {
        winsize.meth = 0
    }

    min.disp=1e-6
    min.InvDisp = 0.001
    max.InvDisp = max(1/max(min.disp, 1e-6), min.InvDisp)

    minMu = 0
    maxMu = 1

    # Get the group labels, THIS ASSUMES CORRECT REFERENCE LEVEL SET
    group2 = levels(pData(meth)[, comparison])[2]
    group1 = levels(pData(meth)[, comparison])[1]

    # Determine which rows of pData belong to which group
    # / which columns of Cov and M matrices belong to which group
    group2_idx = which(pData(meth)[,comparison] == group2)
    group1_idx = which(pData(meth)[,comparison] == group1)

    # Determine which sample column indexes to use for dispersion calculation
    if(dispersion == 'both') {
        disp_idx = c(group2_idx, group1_idx)
    } else if (dispersion == group2) {
        disp_idx = group2_idx
    } else if (dispersion == group1) {
        disp_idx = group1_idx
    } else {
        stop('"dispersion" should be one of "both", the name of group2, or the name of group1')
    }

    all_cov = bsseq::getCoverage(meth, type = 'Cov')
    all_m = bsseq::getCoverage(meth, type = 'M')

    muEst = matrix(0, ncol = ncol(meth), nrow = nrow(meth))
    muEst[, group2_idx] = DelayedArray::rowSums(all_m[, group2_idx]) / (DelayedArray::rowSums(all_cov[, group2_idx]) + 1e-100)
    muEst[, group1_idx] = DelayedArray::rowSums(all_m[, group1_idx]) / (DelayedArray::rowSums(all_cov[, group1_idx]) + 1e-100)
    muEst = DelayedArray::DelayedArray(muEst)

    # Determine which loci satisfy min.per.group
    valid_idx = which(
        DelayedArray::rowSums(bsseq::getCoverage(meth, type = 'Cov')[, group2_idx] > 0) >= min.per.group[2] & DelayedArray::rowSums(bsseq::getCoverage(meth, type = 'Cov')[, group1_idx] > 0) >= min.per.group[1]
    )

    # Go through each index for a valid locus
    results = GRangesList(mclapply(valid_idx[1:10], function(idx){

        # Pull out the BSseq object for just the locus of interest
        locus = meth[idx]

        # Setup GRanges where we expand from each CpG by the winsize
        # Since extending in both directions, divide the winsize.* by 2
        # This is a GRanges object with only 1 range
        disp_win_gr = GenomicRanges::flank(GenomicRanges::granges(locus), width = winsize.disp, start = TRUE, both = TRUE)
        meth_win_gr = GenomicRanges::flank(GenomicRanges::granges(locus), width = winsize.meth, start = TRUE, both = TRUE)

        # If not using local information, then adjust the result of GenomicRanges::flank(). I think
        # there is a bug where if the width is set to 0, it actually subtracts 1 from then end.
        # NOTE: Consider reporting this as a bug to Bioc
        if(winsize.disp == 0){
            end(disp_win_gr) = end(disp_win_gr) + 1
        }
        if(winsize.meth == 0){
            end(meth_win_gr) = end(meth_win_gr) + 1
        }

        # Look for overlaps in the original meth for each of the above two
        # NOTE: You can use meth (BSseq) and it will access the GRanges for GenomicRanges::findOverlaps()
        disp_overlaps = GenomicRanges::findOverlaps(query = meth, subject = disp_win_gr, ignore.strand = TRUE)
        meth_overlaps = GenomicRanges::findOverlaps(query = meth, subject = meth_win_gr, ignore.strand = TRUE)

        # Normalize the locations of the loci in the window to be in the domain of
        # the weight function [-1, 1]
        # locus is GRanges with a single range
        # disp_loci and meth_loci are GRanges with number of ranges equal to overlaps
        # based on the window size
        disp_loci = meth[S4Vectors::queryHits(disp_overlaps)]
        meth_loci = meth[S4Vectors::queryHits(meth_overlaps)]

        # Each is a vector of input values to the weight function
        # We need to scale the loci in the window onto the interval [-1, 1] because
        # that is the domain of the weightFuncction.
        disp_loci_norm = (start(disp_loci) - start(locus)) / width(disp_win_gr)
        meth_loci_norm = (start(meth_loci) - start(locus)) / width(meth_win_gr)

        # Calculate the weights
        # Each is a vector of values of the weight function.
        weightPhi = weightFunc(disp_loci_norm)
        weightMu = weightFunc(meth_loci_norm)

        # Collect Cov and M matrices for all the loci in the window
        # These are delayed matrices. Rows are loci and columns are samples
        locus_cov = bsseq::getCoverage(locus, type = 'Cov')
        locus_meth = bsseq::getCoverage(locus, type = 'M')

        disp_cov = bsseq::getCoverage(meth[S4Vectors::queryHits(disp_overlaps), ], type = 'Cov')
        disp_meth = bsseq::getCoverage(meth[S4Vectors::queryHits(disp_overlaps), ], type = 'M')

        meth_cov = bsseq::getCoverage(meth[S4Vectors::queryHits(meth_overlaps), ], type = 'Cov')
        meth_meth = bsseq::getCoverage(meth[S4Vectors::queryHits(meth_overlaps), ], type = 'M')

        # Collect the correct rows of muEst
        muEst_local = muEst[S4Vectors::queryHits(disp_overlaps),]

        # Convert to old methylSig notion of C reads (methylated) and T reads (unmethylated)
        # Then we can reuse the derivative and log likelihood functions Yongseok implemented.
        # These are delayed matrices. Rows are loci and columns are samples
        disp_local_creads = disp_meth
        disp_local_treads = disp_cov - disp_meth

        meth_local_creads = meth_meth
        meth_local_treads = meth_cov - meth_meth

        # Compute the degrees of freedom for the locus
        if(dispersion == 'both') {
            df_subtract = 2
        } else {
            df_subtract = 1
        }
        df = pmax(DelayedArray::rowSums(locus_cov[, disp_idx] > 0) - df_subtract, 0)
        # Compute the degrees of freedom to be used in the test for differential methylation
        df = sum(df * weightPhi)

        if(df > 1) {
            ### Common dispersion calculation
            # This returns a singleton numeric
            if(methylSig_derivativePhi(
                phi = max.InvDisp,
                lCreads = disp_local_creads[, disp_idx],
                lTreads = disp_local_treads[, disp_idx],
                mu = muEst_local[, disp_idx],
                weight = weightPhi) >= 0) {

                # Describe
                phiCommonEst = max.InvDisp
            } else if(methylSig_derivativePhi(
                phi = min.InvDisp,
                lCreads = disp_local_creads[, disp_idx],
                lTreads = disp_local_treads[, disp_idx],
                mu = muEst_local[, disp_idx],
                weight = weightPhi) <= 0){

                # Describe
                phiCommonEst = min.InvDisp
            } else {
                # Describe
                phiCommonEst = stats::uniroot(methylSig_derivativePhi,
                    c(min.InvDisp, max.InvDisp),
                    disp_local_creads[, disp_idx],
                    disp_local_treads[, disp_idx],
                    muEst_local[, disp_idx],
                    weightPhi)$root
            }

            ### Common group means calculation
            # This returns a numeric vector (group1 and then group2) with the mu estimate
            muEstC = sapply(list(group1_idx, group2_idx, c(group1_idx, group2_idx)), function(group_idx){
                if(sum(meth_local_creads[, group_idx]) == 0) {
                    # If there are no local C reads, methylation is 0
                    return(0)
                } else if (sum(meth_local_treads[, group_idx]) == 0) {
                    # If there are no local T reads, methylation is 1
                    return(1)
                } else {
                    # Otherwise, do something fancier
                    return(
                        stats::uniroot(methylSig_derivativeMu ,
                            c(minMu, maxMu),
                            meth_local_creads[, group_idx],
                            meth_local_treads[, group_idx],
                            phiCommonEst,
                            weightMu)$root
                    )
                }
            })

            ### log Likelihood ratio calculation
            logLikRatio =
                methylSig_logLik(
                    mu = muEstC[1],
                    phi = phiCommonEst,
                    lCreads = meth_local_creads[, group1_idx],
                    lTreads = meth_local_treads[, group1_idx],
                    weight = weightMu) +
                methylSig_logLik(
                    mu = muEstC[2],
                    phi = phiCommonEst,
                    lCreads = meth_local_creads[, group2_idx],
                    lTreads = meth_local_treads[, group2_idx],
                    weight = weightMu) -
                methylSig_logLik(
                    mu = muEstC[3],
                    phi = phiCommonEst,
                    lCreads = meth_local_creads[, c(group1_idx, group2_idx)],
                    lTreads = meth_local_treads[, c(group1_idx, group2_idx)],
                    weight = weightMu)

            locus = GenomicRanges::granges(locus)
            ### Add results to mcols(locus)
            mcols(locus)$phiCommonEst = phiCommonEst
            mcols(locus)$logLikRatio = logLikRatio
            mcols(locus)$muEstC_group1 = muEstC[1]*100
            mcols(locus)$muEstC_group2 = muEstC[2]*100
            mcols(locus)$muEstC_group12 = muEstC[3]*100
            mcols(locus)$df = df + 2
        } else {
            # Not enough degrees of freedom, return NAs
            mcols(locus)$phiCommonEst = NA
            mcols(locus)$logLikRatio = NA
            mcols(locus)$muEstC_group1 = NA
            mcols(locus)$muEstC_group2 = NA
            mcols(locus)$muEstC_group12 = NA
            mcols(locus)$df = NA
        }

        return(locus)
    }, mc.cores = num.cores))

    results = unlist(results)

    if(T.approx) {
         results$pvalue = stats::pt(-sqrt(pmax(results$logLikRatio, 0)), results$df) * 2
    } else {
         results$pvalue = stats::pchisq(pmax(results$logLikRatio, 0), 1, lower.tail = F)
    }

    # Set any methylation difference less than 0.01 to 0
    results$meth.diff = (results$muEstC_group2 - results$muEstC_group1)
    results$meth.diff[abs(results$meth.diff) < 0.01] = 0
    results$meth.diff = as.numeric(results$meth.diff)

    return(results)
}





# Not called by ny other function
methylSigDf <- function(meth, groups=c("Treatment"=1,"Control"=0), min.per.group=c(3,3)) {
    treatment = slot(meth, "treatment")

    group1 = which(treatment == groups[1])
    group2 = which(treatment == groups[2])

    if(length(min.per.group) == 1) {
        min.per.group = c(min.per.group,min.per.group)
    }

    nSamples1 = rowSums(do.call(cbind,meth@.Data[slot(meth,"coverage.index")[group1]]) > 0, na.rm = TRUE)
    nSamples2 = rowSums(do.call(cbind,meth@.Data[slot(meth,"coverage.index")[group2]]) > 0, na.rm = TRUE)

    whichValidList = (nSamples1 >= min.per.group[1] &
                      nSamples2 >= min.per.group[2])

    (nSamples1+nSamples2- 2)[whichValidList]
}

#' Weighted auto correlation of methylation rates
#'
#' This funciton calculates the weighted auto correlation of methylation rates based on the coverage level at each CpG site.
#'
#' The weight for the locus i and j is \code{coverage[i]*coverage[j]/(coverage[i]+coverage[j])}.
#'
#' @param chr A single string indicating chromosome from \code{methylSigData-class} or \code{methylSigDiff-class} objects.
#' @param start A vector of CpG site loci from \code{methylSigData-class} or \code{methylSigDiff-class} objects.
#' @param methRates A numeric vector of methylation rates at each CpG site.
#' @param coverage A numeric vector of coverage levels (number of reads) at each CpG site.
#' @param lags Autocorrelation lags. This argument can be a vector if you want to calculate multiiple lags. Default is 2.
#'
#' @return Auto correlation value or vector for the lags entered.
#'
#' @examples
#' data(sampleData)
#'
#' lags = c(2,5,10,50)
#'
#' autoCorr = methylSigWeightedAutoCorr(meth[,"chr"], meth[,"start"],
#'              meth[,"numCs1"]/meth[,"coverage1"], meth[,"coverage1"],
#'              lags=lags)
#' autoCorr
#'
#' @export
# Not called by any other function
methylSigWeightedAutoCorr<-function(chr, start, methRates, coverage, lags=2) {
    ret = rep(0, NROW(lags))

    MAXSTART = max(start) + max(lags) + 1
    cpgSites = as.numeric(chr)*MAXSTART + start

    ord = order(cpgSites)
    cpgSites = cpgSites[ord]
    methRates = methRates[ord]
    coverage = coverage[ord]

    validList = which(coverage > 0)

    cpgValid = cpgSites[validList]
    methRates = methRates[validList]
    coverage = coverage [validList]

    for(i in 1:length(lags)) {
        lag = lags[i]

        cpgToTest = cpgValid + lag

        cpgNextIndex = findInterval(cpgToTest, cpgValid)
        foundList = which(cpgValid[cpgNextIndex] == cpgToTest)

        ret[i] = corr(cbind(methRates[foundList], methRates[cpgNextIndex[foundList]]),
              w=coverage[foundList]*coverage[cpgNextIndex[foundList]]/(coverage[foundList]+ coverage[cpgNextIndex[foundList]]))
    }

    ret
}
