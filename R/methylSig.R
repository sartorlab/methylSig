#' Write a methylSigDiff object as a text file
#'
#' This funciton writes \code{methylSigDiff-class} as a text file. All options for write.table are avaiable for this function.
#'
#' @param object a \code{methylSigDiff-class} object.
#' @param ... Arguments inherited from \code{\link[utils]{write.table}}
#'
#' @return None
#'
#' @examples
#' data(sampleData)
#' write.methylSigDiff(myDiffSigboth, file="myResult.txt", row.names=FALSE, quote=FALSE, sep="\t")
#'
#' @export
write.methylSigDiff <- function(object, ...) {
    printRange = 1:NROW(object@data.ids)
    printData = data.frame(chr=object@data.chr[printRange], start=object@data.start[printRange],
                           end=object@data.end[printRange], strand=object@data.strand[printRange],
                           object@results[printRange,,drop=FALSE])
    write.table(printData, ...)
}

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

######## Functions ########
#                         #
#    Row: Samples         #
#    Column: Locations    #
###########################

# Called by methylSig_dataProcess
methylSig_derivativePhi <- function(phi, lCreads, lTreads, mu, weight) {
    derivative <- 0
    if(NCOL(lCreads) == 1) {
        ### Only one location, weight does not matter
            vlist <- which(lCreads > 0)
            if(length(vlist) > 0)
                derivative = derivative + sum( mu[vlist] * (digamma((mu[vlist] * phi) + lCreads[vlist]) - digamma(mu[vlist] * phi + 1e-100)) )

            vlist <- which(lTreads > 0)
            if(length(vlist) > 0)
                derivative = derivative + sum( (1 - mu[vlist]) * (digamma( ((1-mu[vlist]) * phi) + lTreads[vlist]) - digamma( ((1-mu[vlist]) * phi) + 1e-100)))

            vlist <- which((lCreads+lTreads) > 0)
            if(length(vlist) > 0)
                derivative = derivative - sum( digamma(phi + lCreads[vlist] + lTreads[vlist]) - digamma(phi))
    } else {
        for(g in 1:NROW(lCreads)) {
            vlist <- which(lCreads[g,] > 0)
            if(length(vlist) > 0)
                derivative = derivative + sum( weight[vlist] * mu[g,vlist] * (digamma(mu[g,vlist] * phi + lCreads[g,vlist]) - digamma(mu[g,vlist] * phi + 1e-100)) )

            vlist <- which(lTreads[g,] > 0)
            if(length(vlist) > 0)
                derivative = derivative + sum( weight[vlist] * (1 - mu[g,vlist]) * (digamma((1 - mu[g,vlist]) * phi + lTreads[g,vlist]) - digamma((1 - mu[g,vlist]) * phi + 1e-100)) )

            vlist <- which((lCreads[g,]+lTreads[g,]) > 0)
            if(length(vlist) > 0)
                derivative = derivative - sum( weight[vlist] * (digamma(phi + lCreads[g,vlist] + lTreads[g,vlist]) - digamma(phi)) )
        }
    }

    derivative
}

# Called by methylSig_dataProcess
methylSig_derivativeMu <- function(mu, lCreads, lTreads, phi, weight) {
    derivative <- 0
    if(NCOL(lCreads) == 1) {
        vlist <- which(lCreads > 0)
        if(length(vlist) > 0)
            derivative = derivative + sum(digamma(mu*phi+lCreads[vlist])-digamma(mu*phi+1e-100))
        vlist <- which(lTreads > 0)
        if(length(vlist) > 0)
            derivative = derivative - sum(digamma((1-mu)*phi+lTreads[vlist])-digamma((1-mu)*phi+1e-100))
    } else {
        for(g in 1:NROW(lCreads)) {
            vlist <- which(lCreads[g,] > 0)
            if(length(vlist) > 0)
                derivative = derivative + sum(weight[vlist]*(digamma(mu*phi+lCreads[g,vlist])-digamma(mu*phi+1e-100)))
            vlist <- which(lTreads[g,] > 0)
            if(length(vlist) > 0)
                derivative = derivative - sum(weight[vlist]*(digamma((1-mu)*phi+lTreads[g,vlist])-digamma((1-mu)*phi+1e-100)))
        }
    }
    derivative
}

# Called by methylSig_dataProcess
methylSig_logLik  <- function(mu, phi, lCreads, lTreads, weight) {
    llik  <- 0
    if(NCOL(lCreads) == 1) {
        ###### single location, weight should be 1
        for(g in 1:NROW(lCreads)) {
            if(lCreads[g] > 0) llik = llik + lgamma(mu*phi+lCreads[g]) - lgamma(mu*phi+1e-100)
            if(lTreads[g] > 0) llik = llik + lgamma((1-mu)*phi+lTreads[g])-lgamma((1-mu)*phi+1e-100)
        }
    } else {
        for(g in 1:NROW(lCreads)) {
            vlist <- which(lCreads[g,] > 0)
            if(length(vlist) > 0) for(i in vlist) { llik = llik + weight[i]*(lgamma(mu*phi+lCreads[g,i])-lgamma(mu*phi+1e-100))}
            vlist <- which(lTreads[g,] > 0)
            if(length(vlist) > 0) for(i in vlist) { llik = llik + weight[i]*(lgamma((1-mu)*phi+lTreads[g,i])-lgamma((1-mu)+1e-100))}
        }
    }

    2*llik
}

##### Weight function: tri-weight
# Called in methylSigCalc
methylSig_weightFunc <- function(u) (1-u^2)^3

# Called by methylSigCalc
# loc is an index for a locus in the obj
methylSig_dataProcess <- function(loc,obj){
    minMu = 0
    maxMu = 1
    group1=obj$groups[[1]]
    group2=obj$groups[[2]]

    ### all Groups is used to calculate common dispersion prameters)
    # NOTE: This is not used
    allGroupsIndex = 3

    # This is really just nLoci again...
    locSize = NCOL(obj$creads)
    ############################

    # Counter
    if(obj$whichOrd[loc] %%10000 == 0) {
       if(obj$whichOrd[loc] %% 50000 == 0) cat("|") else cat(".")
       if(obj$whichOrd[loc] %% 200000 == 0) {
           if(obj$whichOrd[loc] >= 1000000) {
               cat("(",obj$whichOrd[loc] / 1000000, "m)\n", sep="")
           } else {
               cat("(",obj$whichOrd[loc] / 1000, "k)\n", sep="")
           }
       }
    }

    # The wMethIndex and wDispIndex is based on the window
    # Let's say loc = 2, then saying we have 1000 loci:
    # validMuList = max(1, 2 - 100):min(1000, 2 + 100) = 1:102
    validMuList  <- max(1, loc - obj$wMethIndex):min(locSize, loc + obj$wMethIndex )
    validPhiList <- max(1, loc - obj$wDispIndex):min(locSize, loc + obj$wDispIndex )

    # Converts an index into actual locations by determining which are within the window?
    validMuList  = validMuList[which(abs(obj$uniqueLoc[validMuList]  - obj$uniqueLoc[loc]) <= obj$wMeth)]
    validPhiList = validPhiList[which(abs(obj$uniqueLoc[validPhiList] - obj$uniqueLoc[loc]) <= obj$wDispersion)]

    # Of the locations in the window, which can be used for the phi calculation?
    # These are the locations in the windows whose coverage is > 0 (depending on both or either group as the choice for dispersion estimate)
    whichUseful = which(obj$validForPhiCalculate[validPhiList] > 0)
    if(length(whichUseful) == 0)  return(c(loc,NA,NA,NA,NA,NA,NA))

    # Subset the validPhiList by whichUseful
    validPhiList = validPhiList[whichUseful]

    # If there are more than 5 CpGs in validPhiList, restrict to the first 5?
    # These 5 will be the closest ones?
    if(length(validPhiList) > 5) validPhiList = validPhiList[order(abs(obj$uniqueLoc[validPhiList]  - obj$uniqueLoc[loc]))[1:5]]

    # These are the weights
    # NOTE: I don't quite understand the calculation... What is the second half of the difference?
    # The function has domain [-1, 1] so you need to see how far away each thing in the window
    # is from the CpG and then normalize it to be in that domain.
    weightPhi <- obj$weightFunc((obj$uniqueLoc[validPhiList] - obj$uniqueLoc[loc])*obj$wDispNorm)
    weightMu  <- obj$weightFunc((obj$uniqueLoc[validMuList]  - obj$uniqueLoc[loc])*obj$wMethNorm)

#    weightMu = weightMu / sum(weightMu)

    df = sum(obj$validForPhiCalculate[validPhiList]*weightPhi)
    if(df > 1) {
##          && sum(obj$creads[group1,validMuList] +  obj$treads[group1,validMuList] > 0) >= obj$numValidMu[1]
##          && sum(obj$creads[group2,validMuList] +  obj$treads[group2,validMuList] > 0) >= obj$numValidMu[2]) {

        ##### common dispersion calculation
        if(methylSig_derivativePhi(
            phi = obj$max.InvDisp,
            lCreads = obj$creads[obj$dispersionGroups, validPhiList],
            lTreads = obj$treads[obj$dispersionGroups, validPhiList],
            mu = obj$muEst[obj$dispersionGroups, validPhiList],
            weight = weightPhi) >= 0) {

            phiCommonEst = obj$max.InvDisp
        } else if(methylSig_derivativePhi(
            phi = obj$min.InvDisp,
            lCreads = obj$creads[obj$dispersionGroups, validPhiList],
            lTreads = obj$treads[obj$dispersionGroups, validPhiList],
            mu = obj$muEst[obj$dispersionGroups, validPhiList],
            weight = weightPhi) <= 0){

            phiCommonEst = obj$min.InvDisp
        } else {
            phiCommonEst = uniroot(
                methylSig_derivativePhi,
                    c(obj$min.InvDisp, obj$max.InvDisp),
                    obj$creads[obj$dispersionGroups, validPhiList],
                    obj$treads[obj$dispersionGroups, validPhiList],
                    obj$muEst[obj$dispersionGroups, validPhiList],
                    weightPhi)$root
        }

        ##### common group means calculation

        muEstC = rep(0,NROW(obj$groups))
        for(groupsIndex in 1:NROW(obj$groups)) {
           if(sum(obj$creads[obj$groups[[groupsIndex]], validMuList])==0) {
                muEstC[groupsIndex] = 0
           } else if(sum(obj$treads[obj$groups[[groupsIndex]], validMuList])==0) {
                muEstC[groupsIndex] = 1
           } else {
                muEstC[groupsIndex] = uniroot(methylSig_derivativeMu ,c(minMu,maxMu), obj$creads[obj$groups[[groupsIndex]], validMuList],
                                    obj$treads[obj$groups[[groupsIndex]], validMuList],phiCommonEst, weightMu)$root
           }
        }

        #### log Likelihood ratio calculation
        ###### testing ###########
        validMuList = loc
        weightMu = 1
        ##########################

        logLikRatio = methylSig_logLik (muEstC[1], phiCommonEst, obj$creads[obj$groups[[1]], validMuList], obj$treads[obj$groups[[1]], validMuList], weightMu) +
                           methylSig_logLik (muEstC[2], phiCommonEst, obj$creads[obj$groups[[2]], validMuList], obj$treads[obj$groups[[2]], validMuList], weightMu) -
                           methylSig_logLik (muEstC[3], phiCommonEst, obj$creads[obj$groups[[3]], validMuList], obj$treads[obj$groups[[3]], validMuList], weightMu)

        return(c(loc,phiCommonEst,logLikRatio,muEstC,df+2))
    }

    return(c(loc,NA,NA,NA,NA,NA,NA))
}

#' Calculates differential methylation statistics using a Beta-binomial approach.
#'
#' The function calculates differential methylation statistics between two groups of samples. This is the main function of the methylSig package, and the method most users should use to test for DMCs or DMRs. The function uses a Beta-binomial approach to calculate differential methylation statistics, accounting for variation among samples within each group.
#'
#' The function calculates differential methylation statistics between two groups of samples. The function uses Beta-binomial approach to calculate differential methylation statistics, accounting for variation among samples within each group. Users who wish to tile their data and test for differentially methylated regions (DMRs) instead DMCs should first use the \code{\link{methylSigTile}} function before using this function.
#'
#' @param meth A \code{methylSigData-class} object to calculate differential methylation statistics. It can be obtained using `methylSigReadData'.
#' @param groups A vector of two numbers specify two groups to compare. See `treatment' argument of \code{\link{methylSigReadData}} function. Default is \code{c(Treatment=1,Control=0)}.
#' @param dispersion A value indicating which group or groups are used to estimate the variance (dispersion). If groups are defined as c(Treatment=1,Control=0), dispersion can take values "Treatment", "Control", 1, 0 or "both". Default is "both".
#' @param local.disp A logical value indicating whether to use local information to improve dispersion parameter estimation. Default is FALSE.
#' @param winsize.disp A number to specify the window size in basepairs for local dispersion estimation. The dispersion (variance) parameter of the groups at the particular location LOC is calculated using information from LOC - winsize.disp to LOC + winsize.disp. This argument is only activated when local.disp=TRUE. Default is 200.
#' @param local.meth A logical value indicating whether to use local information to improve methylation level estimation. Default is FALSE.
#' @param winsize.meth a number to specify the window size in basepairs for local methylation level estimation. The group methylation level at the particular location LOC is calculated using information from LOC - winsize.meth to LOC + winsize.meth. This argument is only activated when local.meth=TRUE. Default is 200.
#' @param min.per.group A vector with two numbers that specify the minimum numbers of samples required to be qualify as defferentially methylated region.  If it is a single number, both groups will use it as the minimum requried number of samples. Default is c(3,3).
#' @param weightFunc A weight kernel function. The input of this function is from -1 to 1. The default is the tri-weight kernel function defined as function(u) = (1-u^2)^3. Function value and range of parameter for weight function should be from 0 to 1.
#' @param T.approx A logical value indicating whether to use squared t approximation for the likelihood ratio statistics. Chi-square approximation (T.approx = FALSE) is recommended when the sample size is large.  Default is TRUE.
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
methylSigCalc = function(meth, groups=c("Treatment"=1,"Control"=0), dispersion="both",
         local.disp=FALSE, winsize.disp=200,
         local.meth=FALSE, winsize.meth=200,
         min.per.group=c(3,3), weightFunc=methylSig_weightFunc, T.approx = TRUE,
         num.cores = 1) {

    # Smallest dispersion?
    min.disp=1e-6

    if(meth@resolution == "tfbs") {
        local.disp=FALSE
        local.meth=FALSE
    }

    # Set window sizes to 0 if either is FALSE
    if(local.meth == FALSE) winsize.meth = 0
    if(local.disp == FALSE) winsize.disp = 0

    # Get the treatment vector from the meth object
    treatment = slot(meth,"treatment")

    # Determine the indices in the treatment vector belonging to each group
    group1 = which(treatment == groups[1])
    group2 = which(treatment == groups[2])

    if(length(group1) == 0 || length(group2) == 0) {
        stop("Groups do not match your treatment in the input data")
    }

    # If min.per.group is a single number, make it apply to both groups
    if(length(min.per.group) == 1) {
        min.per.group = c(min.per.group,min.per.group)
    }

    # Instantiate a new environment in which to do the DM testing...?
    # globalenv() is probably not recommended...
    methSigObject=new.env(parent=globalenv())
    class(methSigObject)='pointer'

    #
    methSigObject$groups <- list(group1 = group1,
               group2 = group2,
               group3 = c(group1,group2)
              )

    # Order the data (data.ids takes the form chr1.390, for example)
    orderMethStart = order(meth@data.ids)

    ########################################################
    # Take the transpose of the numTs and numTs from the meth object.
    # Why take the transpose?
    # NOTE: These are indeed counts
    methSigObject$treads   = t(meth@data.numTs[orderMethStart,])
    methSigObject$creads   = t(meth@data.numCs[orderMethStart,])

    # Reset any NAs to 0
    methSigObject$treads[is.na(methSigObject$treads)] = 0
    methSigObject$creads[is.na(methSigObject$creads)] = 0

    ########################################################
    # vector: uniqueLoc is like a hash of the chromosome locations (they are integers)
    methSigObject$uniqueLoc = meth@data.ids[orderMethStart]

    ########################################################
    # integer: Not really sure what to make of this...
    methSigObject$stepSize = ifelse(meth@destranded, 2, 1)

    ############################
    # integer: Window size for for localalized methylation estimates
    methSigObject$wMeth = winsize.meth
    # numeric: If destranded, then shrink the window by half because you're
    # incorporating twice the amount of information?
    # So if winsize.meth = 200 and destranded, then wMethIndex = 100
    methSigObject$wMethIndex = winsize.meth/methSigObject$stepSize
    # numeric: A normalization factor of some sort...?
    # Continuing the above, we'd have 1 / (2 * (100 + 1)) = 0.00495
    methSigObject$wMethNorm = 1/(methSigObject$stepSize*(methSigObject$wMethIndex+1))

    ############################
    # integer: Similar procedure to three preceding lines
    methSigObject$wDispersion = winsize.disp
    # numeric: If destranded, then shrink the window by half because you're
    # incorporating twice the amount of information?
    methSigObject$wDispIndex = winsize.disp/methSigObject$stepSize
    # numeric: A normalization factor of some sort...?
    methSigObject$wDispNorm  = 1/(methSigObject$stepSize*(methSigObject$wDispIndex+1))

    ########################################################
    # Inverted dispersion?
    # numeric:
    methSigObject$min.InvDisp = 0.001
    # numeric: Will always be 1e+6..., so ... why?
    methSigObject$max.InvDisp = max(1/max(min.disp,1e-6), methSigObject$min.InvDisp)

    ########################################################
    methSigObject$numValidMu = min.per.group
    methSigObject$weightFunc = weightFunc

    if(dispersion == "both") {
        methSigObject$dispersionGroups = c(group1,group2)
        # This is a vector of which CpGs will have phi (dispersion) calcualted for them.
        # NOTE: creads and treads are count matrices
        # numCs and numTs are misnomers until fixed in methylSigReadData (old version)
        # df$numCs = round(df$numCs * df$coverage / 100)
        # df$numTs = round(df$numTs * df$coverage / 100)
        # NOTE: Not sure why 1 is subtracted from the colSums...
        # NOTE: Columns are locations and rows are samples
        # locations valid for phi calculation are those with coverage in the both groups (or each, depending) is greater than 0.
        methSigObject$validForPhiCalculate = pmax(colSums((methSigObject$creads + methSigObject$treads > 0)[group2,]) - 1, 0) + pmax(colSums((methSigObject$creads + methSigObject$treads > 0)[group1,]) - 1, 0)
    }  else if(dispersion == groups[1] || (length(names(groups)[1])>0 && dispersion == names(groups)[1])) {
        methSigObject$dispersionGroups = group1
        # Use only group1
        methSigObject$validForPhiCalculate = pmax(colSums((methSigObject$creads + methSigObject$treads > 0)[group1,]) - 1, 0)
    } else if(dispersion == names(groups)[2] || dispersion == groups[2]) {
        methSigObject$dispersionGroups = group2
        # Use only group2
        methSigObject$validForPhiCalculate = pmax(colSums((methSigObject$creads + methSigObject$treads > 0)[group2,]) - 1, 0)
    } else {
        cat("Dispersion should be \"", names(groups)[1], "\", \"", names(groups)[2], "\" or \"both\".\n", sep="")
        return(NULL)
    }

    # For some reason Yongseok took the transpose, so columns are sites and rows are samples...
    nLoci = NCOL(methSigObject$creads)
    # Instantiate a matrix to hold the mu estimates (methylation rates)
    methSigObject$muEst <- matrix(0, ncol=nLoci, nrow=NROW(methSigObject$creads))

    # Get the methylation rates at each site by taking creads / coverages summed over the samples in each group
    muList1 <- colSums(methSigObject$creads[group1,])/(colSums((methSigObject$creads+methSigObject$treads)[group1,])+1e-100)
    muList2 <- colSums(methSigObject$creads[group2,])/(colSums((methSigObject$creads+methSigObject$treads)[group2,])+1e-100)

    # Populate the muEst
    for(g in group1) {
        methSigObject$muEst[g,] = muList1
    }
    for(g in group2) {
        methSigObject$muEst[g,] = muList2
    }

    # Which sites meet the min.per.group threshold?
    validLoci = ((colSums((methSigObject$creads + methSigObject$treads> 0)[group1,]) >= min.per.group[1])
               & (colSums((methSigObject$creads + methSigObject$treads> 0)[group2,]) >= min.per.group[2]))

    # Used for status updates.
    methSigObject$whichOrd = cumsum(validLoci)

    # Some reporting
    nLoci = sum(validLoci)
    if(nLoci >= 1000000) {
        cat("Total number of ", meth@resolution, "s: ", round(nLoci/1000000,2), "m\n", sep="")
    } else if(nLoci >= 1000) {
        cat("Total number of ", meth@resolution, "s: ", round(nLoci/1000,2), "k\n", sep="")
    } else {
        cat("Total number of ", meth@resolution, "s: ", nLoci, "\n", sep="")
    }

    # Go through each validLoci index (loc) and perform methylSig_dataProcess with the methSigObject (obj)
    if(num.cores == 1) {
        result = do.call(rbind, lapply(which(validLoci), methylSig_dataProcess, methSigObject))
    } else {
        result = do.call(rbind, mclapply(which(validLoci), methylSig_dataProcess, methSigObject, mc.cores=num.cores))
    }
    cat("\n")

    logLikRatio = result[,3]

    if(T.approx) {
         pvalue = pt(-sqrt(pmax(logLikRatio,0)),result[,7])*2
    } else {
         pvalue = pchisq(pmax(logLikRatio,0), 1, lower.tail=F)
    }

    # Set any methylation difference less than 0.01 to 0
    meth.diff = (result[,4] - result[,5])*100
    meth.diff[which(abs(meth.diff) < 0.01)] = 0
    meth.diff = as.numeric(meth.diff)

    results=cbind(pvalue,p.adjust(pvalue, method ="BH"), meth.diff, logLikRatio,
                  result[,2], result[,7], result[,4]*100, result[,5]*100)

    colnames(results) = c("pvalue","qvalue", "meth.diff","logLikRatio","theta", "df", paste("mu", groups, sep=""))


    optionForLocalDisp = ifelse(local.disp, paste(" & winsize.disp=", winsize.disp, sep=""), "")
    optionForLocalMeth = ifelse(local.meth, paste(" & winsize.meth=", winsize.meth, sep=""), "")

    options = paste("dispersion=", dispersion, " & local.disp=", local.disp, optionForLocalDisp,
                             " & local.meth=", local.meth, optionForLocalMeth, "& min.per.group=c(",min.per.group[1],
                             ",",min.per.group[2], ")& Total: ", NROW(results), sep="")

    methylSig.newDiff(meth@data.ids[orderMethStart[validLoci]], meth@data.chr[orderMethStart[validLoci]],
                              as.integer(meth@data.start[orderMethStart[validLoci]]),as.integer(meth@data.end[orderMethStart[validLoci]]),
                              meth@data.strand[orderMethStart[validLoci]], results, sample.ids=meth@sample.ids[c(group1,group2)],
                              sample.filenames=meth@sample.filenames[c(group1,group2)],
                              treatment=meth@treatment[c(group1,group2)], destranded=meth@destranded,
                              resolution=meth@resolution,
                              options=options, data.options = meth@options)
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
