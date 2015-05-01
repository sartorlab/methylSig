#' Write a methylSigDiff object as a text file
#'
#' This funciton writes \code{\link{methylSigDiff}} as a text file. All options for write.table are avaiable for this function.
#'
#' @param object a \code{link{methylSigDiff}} object.
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
#' @param meth A \code{\link{methylSigData-class}} object to calculate differential methylation statistics. It can be obtained using \code{\link{methylSigReadData}}.
#' @param groups A vector of two numbers specify two groups to compare. See \code{treatment} argument of \code{\link{methylSigReadData}} function. Default is \code{c(Treatment=1,Control=0)}.
#' @param min.per.group A vector with two numbers that specify the minimum numbers of samples required to be qualify as defferentially methylation region. If it is a single number, both groups will use it as the minimum requried number of samples. Default is \code{c(3,3)}.
#'
#' @return A \code{\link{methylSigDiff-class}} object that contains the differential methylation statistics and chromosomal locations. \code{p.adjust} with \code{method="BH"} option is used for P-value correction.
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

#' Obtain tiled methylation data in non-overlapping continuous windows.
#'
#' This funciton tiles data within windows of a given width across genome. It gives total number of methylated cytosines and total number of reads (coverage) within each tiled window. This is used prior to the \code{\link{methylSigCalc}} function when the user prefers to conduct a tiled analysis instead of a base specific analysis for differential methylation. Tiling may provide higher power to detect significant differences, especially for experiments with low coverage.
#'
#' @param meth A \code{\link{methylSigData-class}} object used to tile data.
#' @param win.size An integer value indicating the desired window size in bps. Default is 25.
#'
#' @return A \code{\link{methylSigData-class}} object.
#'
#' @examples
#' data(sampleData)
#' methTile = methylSigTile(meth)
#' @export
methylSigTile <- function(meth,win.size=25) {
    if(meth@resolution == "region") stop("Object has already been tiled")

    MAXBASE = max(meth@data.start) + win.size + 1
#    MAXBASE10 = MAXBASE = 10^{ceiling(log10(MAXBASE + 1))}

    startList = as.numeric(meth@data.chr)*MAXBASE + meth@data.start
    uniqueStartList = unique(startList - (startList - 1) %% win.size)

    whichRow = match(startList - (startList-1) %% win.size, uniqueStartList)
    coverage = numTs = numCs = matrix(0,nrow=NROW(uniqueStartList), ncol=NCOL(meth@data.coverage))

    for(i in 1:NROW(startList)) {
         coverage[whichRow[i],] = coverage[whichRow[i],] + meth@data.coverage[i,]
         numCs[whichRow[i],] = numCs[whichRow[i],] + meth@data.numCs[i,]
         numTs[whichRow[i],] = numTs[whichRow[i],] + meth@data.numTs[i,]
    }

    methylSig.newData(data.ids=uniqueStartList, data.chr=as.factor(levels(meth@data.chr)[as.integer(uniqueStartList/MAXBASE)]),
                      data.start=uniqueStartList%%MAXBASE, data.end=uniqueStartList%%MAXBASE + win.size - 1,
                      data.strand=factor(rep("*",NROW(uniqueStartList))), data.coverage = coverage, data.numTs = numTs, data.numCs = numCs,
                      sample.ids=meth@sample.ids, treatment=meth@treatment, destranded=meth@destranded,
                      resolution="region", sample.filenames=meth@sample.filenames,options=meth@options)
}

#' Obtain tiled methylation data by tiling (pooling) data that a particular transcription factor (TF) is predicted to bind.
#'
#' Gives total number of methylated cytosines and total number of reads (coverage) for each TF. This is used prior to the \code{\link{methylSigCalc}} function when the user prefers to conduct tests to identify significant level of hypermethylation or hypomethylation across their binding sites.
#'
#' @param meth A \code{\link{methylSigData-class}} object used to tile data.
#' @param tfbsInfo An object that contains transcription factor binding sites information.
#'
#' @return A \code{\link{methylSigData-class}} object.
#'
#' @seealso \code{\link{getTFBSInfo}} and \code{\link{methylSigTile}}.
#' @export
methylSigTileTFBS <- function(meth, tfbsInfo) {
    totalTFBS = length(levels(tfbsInfo[[2]]$name))
    numCs = numTs = matrix(0, ncol=NCOL(meth@data.coverage), nrow=totalTFBS)

    for(chr in levels(meth@data.chr)) {
        tfbsChrId = which(tfbsInfo[[1]]$chr == as.character(chr))
        if(length(tfbsChrId) > 0) {
            tfbsIndexList = tfbsInfo[[1]]$start[tfbsChrId]:tfbsInfo[[1]]$end[tfbsChrId]
            endList = tfbsInfo[[2]]$chromEnd[tfbsIndexList]
            startList = tfbsInfo[[2]]$chromStart[tfbsIndexList]

            whichVlist = which(meth@data.chr==chr)
            whichEND = findInterval(endList, meth@data.start[whichVlist])
            whichSTART = findInterval(startList, meth@data.start[whichVlist]+1)

            for(i in which(whichEND > whichSTART)) {
                whichTFBS = as.numeric(tfbsInfo[[2]]$name[tfbsIndexList[i]])
                numCs[whichTFBS,] = numCs[whichTFBS,] + colSums(meth@data.numCs[whichVlist[(whichSTART[i]+1):whichEND[i]],,drop=F], na.rm=T)
                numTs[whichTFBS,] = numTs[whichTFBS,] + colSums(meth@data.numTs[whichVlist[(whichSTART[i]+1):whichEND[i]],,drop=F], na.rm=T)
            }
        }
    }
    methylSig.newData(data.ids=1:totalTFBS, data.chr=as.factor(levels(tfbsInfo[[2]]$name)),
                      data.start=rep(0, totalTFBS), data.end=rep(0, totalTFBS),
                      data.strand=rep(as.factor("NA"), totalTFBS), data.coverage = numTs+numCs, data.numTs = numTs, data.numCs = numCs,
                      sample.ids=meth@sample.ids, treatment=meth@treatment, destranded=meth@destranded,
                      resolution="TF", sample.filenames=meth@sample.filenames,options="")
}


######## Functions ########
#                         #
#    Row: Samples         #
#    Column: Locations    #
###########################

# Called by methylSig_dataProcess
methylSig_derivativePhi <- function(phi,lCreads,lTreads,mu,weight) {
    derivative <- 0
    if(NCOL(lCreads) == 1) {
        ### Only one location, weight does not matter
            vlist <- which(lCreads > 0)
            if(length(vlist) > 0)
                derivative = derivative + sum(mu[vlist]*(digamma(mu[vlist]*phi+lCreads[vlist])-digamma(mu[vlist]*phi+1e-100)))
            vlist <- which(lTreads > 0)
            if(length(vlist) > 0)
                derivative = derivative + sum((1-mu[vlist])*(digamma((1-mu[vlist])*phi+lTreads[vlist])-digamma((1-mu[vlist])*phi+1e-100)))
            vlist <- which((lCreads+lTreads) > 0)
            if(length(vlist) > 0)
                derivative = derivative - sum(digamma(phi+lCreads[vlist]+lTreads[vlist])-digamma(phi))
    } else {
        for(g in 1:NROW(lCreads)) {
            vlist <- which(lCreads[g,] > 0)
            if(length(vlist) > 0)
                derivative = derivative + sum(weight[vlist]*mu[g,vlist]*(digamma(mu[g,vlist]*phi+lCreads[g,vlist])-digamma(mu[g,vlist]*phi+1e-100)))
            vlist <- which(lTreads[g,] > 0)
            if(length(vlist) > 0)
                derivative = derivative + sum(weight[vlist]*(1-mu[g,vlist])*(digamma((1-mu[g,vlist])*phi+lTreads[g,vlist])-digamma((1-mu[g,vlist])*phi+1e-100)))
            vlist <- which((lCreads[g,]+lTreads[g,]) > 0)
            if(length(vlist) > 0)
                derivative = derivative - sum(weight[vlist]*(digamma(phi+lCreads[g,vlist]+lTreads[g,vlist])-digamma(phi)))
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
methylSig_logLik  <- function(mu,phi,lCreads, lTreads, weight) {
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
methylSig_dataProcess <- function(loc,obj){
    minMu = 0
    maxMu = 1
    group1=obj$groups[[1]]
    group2=obj$groups[[2]]

    ### all Groups is used to calculate common dispersion prameters)
    allGroupsIndex = 3

    locSize = NCOL(obj$creads)
    ############################

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

    validMuList  <- max(1,loc-obj$wMethIndex):min(locSize,loc+obj$wMethIndex)
    validPhiList <- max(1,loc-obj$wDispIndex):min(locSize,loc+obj$wDispIndex)

    validMuList  = validMuList [which(abs(obj$uniqueLoc[validMuList]  - obj$uniqueLoc[loc]) <= obj$wMeth)]
    validPhiList = validPhiList[which(abs(obj$uniqueLoc[validPhiList] - obj$uniqueLoc[loc]) <= obj$wDispersion)]

    whichUseful = which(obj$validForPhiCalculate[validPhiList] > 0)
    if(length(whichUseful) == 0)  return(c(loc,NA,NA,NA,NA,NA,NA))

    validPhiList = validPhiList[whichUseful]

    if(length(validPhiList) > 5) validPhiList = validPhiList[order(abs(obj$uniqueLoc[validPhiList]  - obj$uniqueLoc[loc]))[1:5]]

    weightPhi <- obj$weightFunc((obj$uniqueLoc[validPhiList] - obj$uniqueLoc[loc])*obj$wDispNorm)
    weightMu  <- obj$weightFunc((obj$uniqueLoc[validMuList]  - obj$uniqueLoc[loc])*obj$wMethNorm)

#    weightMu = weightMu / sum(weightMu)

    df = sum(obj$validForPhiCalculate[validPhiList]*weightPhi)
    if(df > 1) {
##          && sum(obj$creads[group1,validMuList] +  obj$treads[group1,validMuList] > 0) >= obj$numValidMu[1]
##          && sum(obj$creads[group2,validMuList] +  obj$treads[group2,validMuList] > 0) >= obj$numValidMu[2]) {

        ##### common dispersion calculation
        if(methylSig_derivativePhi(obj$max.InvDisp,obj$creads[obj$dispersionGroups, validPhiList],
                            obj$treads[obj$dispersionGroups,validPhiList],obj$muEst[obj$dispersionGroups,validPhiList], weightPhi) >= 0) {
            phiCommonEst = obj$max.InvDisp
        } else if(methylSig_derivativePhi(obj$min.InvDisp,obj$creads[obj$dispersionGroups,validPhiList],
                            obj$treads[obj$dispersionGroups,validPhiList],obj$muEst[obj$dispersionGroups,validPhiList], weightPhi) <= 0){
            phiCommonEst = obj$min.InvDisp
        } else {
            phiCommonEst = uniroot(methylSig_derivativePhi,c(obj$min.InvDisp,obj$max.InvDisp), obj$creads[obj$dispersionGroups,validPhiList],
                            obj$treads[obj$dispersionGroups,validPhiList],obj$muEst[obj$dispersionGroups,validPhiList], weightPhi)$root
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

        return(c(loc,phiCommonEst,logLikRatio,muEstC,df))
    }

    return(c(loc,NA,NA,NA,NA,NA,NA))
}

#' Calculates differential methylation statistics using a Beta-binomial approach.
#'
#' The function calculates differential methylation statistics between two groups of samples. This is the main function of the methylSig package, and the method most users should use to test for DMCs or DMRs. The function uses a Beta-binomial approach to calculate differential methylation statistics, accounting for variation among samples within each group.
#'
#' The function calculates differential methylation statistics between two groups of samples. The function uses Beta-binomial approach to calculate differential methylation statistics, accounting for variation among samples within each group. Users who wish to tile their data and test for differentially methylated regions (DMRs) instead DMCs should first use the \code{\link{methylSigTile}} function before using this function.
#'
#' @param meth A \code{\link{methylSigData-class}} object to calculate differential methylation statistics. It can be obtained using `methylSigReadData'.
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
#' @return `methylSigDiff' object containing the differential methylation statistics and locations. p.adjust with method="BH" option is used for P-value correction.
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

    ###### future properties
    min.disp=1e-6
    ######

    if(meth@resolution == "tfbs") {
        local.disp=FALSE
        local.meth=FALSE
    }

    if(local.meth == FALSE) winsize.meth = 0
    if(local.disp == FALSE) winsize.disp = 0

    treatment = slot(meth,"treatment")

    group1 = which(treatment == groups[1])
    group2 = which(treatment == groups[2])

    if(length(group1) == 0 || length(group2) == 0) {
        stop("Groups do not match your treatment in the input data")
    }

    if(length(min.per.group) == 1) {
        min.per.group = c(min.per.group,min.per.group)
    }

    methSigObject=new.env(parent=globalenv())
    class(methSigObject)='pointer'
    methSigObject$groups <- list(group1 = group1,
               group2 = group2,
               group3 = c(group1,group2)
              )

    orderMethStart = order(meth@data.ids)

    methSigObject$treads   = t(meth@data.numTs[orderMethStart,])
    methSigObject$creads   = t(meth@data.numCs[orderMethStart,])
    methSigObject$treads[is.na(methSigObject$treads)] = 0
    methSigObject$creads[is.na(methSigObject$creads)] = 0

    methSigObject$uniqueLoc = meth@data.ids[orderMethStart]
    methSigObject$stepSize = ifelse(meth@destranded, 2, 1)
    methSigObject$wMeth = winsize.meth
    methSigObject$wMethIndex = winsize.meth/methSigObject$stepSize
    methSigObject$wMethNorm = 1/(methSigObject$stepSize*(methSigObject$wMethIndex+1))
    methSigObject$wDispersion = winsize.disp
    methSigObject$wDispIndex = winsize.disp/methSigObject$stepSize
    methSigObject$wDispNorm  = 1/(methSigObject$stepSize*(methSigObject$wDispIndex+1))
    methSigObject$min.InvDisp = 0.001
    methSigObject$max.InvDisp = max(1/max(min.disp,1e-6), methSigObject$min.InvDisp)
    methSigObject$numValidMu = min.per.group
    methSigObject$weightFunc = weightFunc

    if(dispersion == "both") {
       methSigObject$dispersionGroups = c(group1,group2)
       methSigObject$validForPhiCalculate = pmax(colSums((methSigObject$creads + methSigObject$treads> 0)[group2,]) - 1, 0) + pmax(colSums((methSigObject$creads + methSigObject$treads> 0)[group1,]) - 1, 0)
    }  else if(dispersion == groups[1] || (length(names(groups)[1])>0 && dispersion == names(groups)[1])) {
        methSigObject$dispersionGroups = group1
        methSigObject$validForPhiCalculate = pmax(colSums((methSigObject$creads + methSigObject$treads > 0)[group1,]) - 1, 0)
    } else if(dispersion == names(groups)[2] || dispersion == groups[2]) {
        methSigObject$dispersionGroups = group2
        methSigObject$validForPhiCalculate = pmax(colSums((methSigObject$creads + methSigObject$treads> 0)[group2,]) - 1, 0)
    } else {
        cat("Dispersion should be \"", names(groups)[1], "\", \"", names(groups)[2], "\" or \"both\".\n", sep="")
        return(NULL)
    }

    nLoci = NCOL(methSigObject$creads)
    methSigObject$muEst <- matrix(0, ncol=nLoci, nrow=NROW(methSigObject$creads))

    muList1 <- colSums(methSigObject$creads[group1,])/(colSums((methSigObject$creads+methSigObject$treads)[group1,])+1e-100)
    muList2 <- colSums(methSigObject$creads[group2,])/(colSums((methSigObject$creads+methSigObject$treads)[group2,])+1e-100)

    for(g in group1) {
        methSigObject$muEst[g,] = muList1
    }
    for(g in group2) {
        methSigObject$muEst[g,] = muList2
    }

    validLoci = ((colSums((methSigObject$creads + methSigObject$treads> 0)[group1,]) >= min.per.group[1])
               & (colSums((methSigObject$creads + methSigObject$treads> 0)[group2,]) >= min.per.group[2]))

    methSigObject$whichOrd = cumsum(validLoci)

    nLoci = sum(validLoci)
    if(nLoci >= 1000000) {
        cat("Total number of ", meth@resolution, "s: ", round(nLoci/1000000,2), "m\n", sep="")
    } else if(nLoci >= 1000) {
        cat("Total number of ", meth@resolution, "s: ", round(nLoci/1000,2), "k\n", sep="")
    } else {
        cat("Total number of ", meth@resolution, "s: ", nLoci, "\n", sep="")
    }

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

    results=cbind(pvalue,p.adjust(pvalue, method ="BH"), (result[,4] - result[,5])*100, logLikRatio,
                  result[,2], result[,7], result[,4]*100, result[,5]*100)

    colnames(results) = c("pvalue","qvalue", "meth.diff","logLikRatio","theta", "df", paste("mu", groups, sep=""))


    optionForLocalDisp = ifelse(local.disp, paste(" & winsize.disp=", winsize.disp, sep=""), "")
    optionForLocalMeth = ifelse(local.meth, paste(" & winsize.meth=", winsize.meth, sep=""), "")

    options = paste("dispersion=", dispersion, " & local.disp=", local.disp, optionForLocalDisp,
                             " & local.meth=", local.meth, optionForLocalMeth, "& min.per.group=c(",min.per.group[1],
                             ",",min.per.group[2], ")& Total: ", NROW(results), sep="")

    methylSig.newDiff(meth@data.ids[orderMethStart[validLoci]], meth@data.chr[orderMethStart[validLoci]],
                              meth@data.start[orderMethStart[validLoci]],meth@data.end[orderMethStart[validLoci]],
                              meth@data.strand[orderMethStart[validLoci]], results, sample.ids=meth@sample.ids[c(group1,group2)],
                              sample.filenames=meth@sample.filenames[c(group1,group2)],
                              treatment=meth@treatment[c(group1,group2)], destranded=meth@destranded,
                              resolution=meth@resolution,
                              options=options, data.options = meth@options)
}

#' Data visualization function
#'
#' Generates data visualization plot of methylation data for a specified genomic interval.
#'
#'   This function offers a unique two-tiered visualization of the methylation data depending on the zoom level. For narrow regions (<1mbp) where at most 500 CpG sites have data reads, users can visualize sample-specific coverage levels and percent methylation at each site, together with group averages, significance levels and a number of genomic annotations.
#'
#' @param meth A `methylSigData' object.
#' @param chr Chromosome in character that matches in methylation data (meth).
#' @param loc.range A vector of two numbers (from, to) to specify the region to visualize on chromosome `chr'.
#' @param groups A vector of two numbers to specify two groups to show in the plot.
#' @param cpgInfo CpG island information to annotate. If missing, no CpG island information will be shown.
#' @param refGeneInfo refGene information to annotate. If missing, no refGene island information will be shown.
#' @param myDiff Differential methylation analysis results from `methylSigCalc'. If missing, no ``-log10[pvalues]'' will be shown.
#' @param tfbsInfo Transcription factor binding site informaiton from `getTFBSInfo' function. If missing, no tfbs will be shown.
#' @param noGroupEst If noGroupEst = TRUE, group estimates will not be drawn in the plot. Default is FALSE
#' @param noDataLine If noDataLine = FALSE, a line will link the same letter that represents the same sample. Default is TRUE.
#' @param cex A numerical value giving the amount by which plotting letters of data should be magnified relative to the default.
#' @param tfbsDense If tfbsDense = TRUE, the plot of transcription factor binding sites is using dense mode. If the region is broader, the tfbsDense is automatically set to TRUE.
#' @param sigQ Maximum significant q-value to be drawn as red line. Default is 0.05.
#'
#' @examples
#' data(sampleData)
#'
#' cpgInfo = getCpGInfo(system.file("annotation", "cpgi.hg18.bed.txt",
#'                  package = "methylSig"))
#' tfbsInfo = getTFBSInfo(system.file("annotation", "tfbsUniform.txt",
#'                  package = "methylSig"))
#' refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt",
#'                  package = "methylSig"))
#'
#' methylSigPlot(meth, "chr21", c(43000000, 44000000), groups=c(1,0),
#'              cpgInfo=cpgInfo, refGeneInfo=refGeneInfo,
#'              myDiff=myDiffSigboth, tfbsInfo=tfbsInfo,
#'              tfbsDense=FALSE, sigQ=0.05)
#'
#' methylSigPlot(meth, "chr21", c(43800000, 43900000), groups=c(1,0),
#'              cpgInfo=cpgInfo, refGeneInfo=refGeneInfo,
#'              myDiff=myDiffSigboth, tfbsInfo=tfbsInfo,
#'              tfbsDense=FALSE, sigQ=0.05)
#'
#' @export
methylSigPlot <-function(meth, chr, loc.range, groups=c("Treatment"=1,"Control"=0), cpgInfo, refGeneInfo, myDiff, tfbsInfo,
             noGroupEst = FALSE, noDataLine = TRUE, cex=1, tfbsDense=TRUE, sigQ=0.05) {
    treatment = slot(meth,"treatment")

    group1 = which(treatment == groups[1])
    group2 = which(treatment == groups[2])
    groupAll = c(group1, group2)

    if(length(group1) == 0 || length(group2) == 0) {
        stop("Groups do not match your treatment in the input data")
    }

    pchList <- c("A", "B", "C", "D", "E", "F", "G", "H", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "2", "3", "4", "5", "6", "7", "8", "9")

    tooWide = 0
    extra = 0

    validList <- which((meth[,"chr"] == chr & meth[,"start"] >= loc.range[1]) & (meth[,"start"] <= loc.range[2]))
    if(diff(loc.range)>= 1000000 || NROW(validList) > 500) {
        warning("Range of the drawing area is two wide, please reduce the range!")
        tfbsDense=TRUE
        noGroupEst = TRUE
        tooWide = 1
    }

    par(mar=c(4,6.5,2,2))

    EXTRA = 12
    EXTRA_TFBS = 30
    GAPS_DATA_ANNOT = 10

    if(!missing(cpgInfo)) extra = extra - EXTRA
    if(!missing(refGeneInfo)) extra = extra - EXTRA
    if(!missing(myDiff)) extra = extra - EXTRA
    if(!missing(tfbsInfo)) extra = ifelse(tfbsDense, extra-EXTRA, extra - EXTRA_TFBS)

    if(tooWide) {
        plot(loc.range,c(0,0), yaxt="n", type="n", xlim = loc.range,ylim=c(extra,0),lwd=3,xlab="",ylab="", col=2)
        yLLoc = -GAPS_DATA_ANNOT
        yHLoc = yLLoc+GAPS_DATA_ANNOT
    } else {
        #extra = extra - GAPS_DATA_ANNOT
        plot(loc.range,c(0,0), yaxt="n", type="n", xlim = loc.range,ylim=c(extra,100),lwd=3,xlab="",ylab="methylation rate", col=2)
        axis(2, at=c(0:5)*20, las=2)
        yLLoc = -EXTRA-GAPS_DATA_ANNOT
        yHLoc = yLLoc+GAPS_DATA_ANNOT
    }

    ###############
    ###############
    ###############
    CPG_ISLAND_COLOR="chartreuse4"
    CPG_SHORE_COLOR="green3"
    CPG_SHELVE_COLOR="gold2"
    CPG_SEA_COLOR="BLUE"
    ###############
    ###############
    ###############

    if(!tooWide && !missing(cpgInfo)) {
        whichChr = which(cpgInfo[[1]]$chr == chr)
        whichRange = cpgInfo[[1]]$start[whichChr]:cpgInfo[[1]]$end[whichChr]
        whichUsed = whichRange[which(cpgInfo[[2]]$start[whichRange] > loc.range[1]-4000 & cpgInfo[[2]]$end[whichRange] < loc.range[2] + 4000)]

        cpgInfoUsed = data.frame(start = c(-Inf, cpgInfo[[2]]$start[whichUsed], Inf), end = c(-Inf, cpgInfo[[2]]$end[whichUsed], Inf))
        cpgId = findInterval(meth[validList,"start"], cpgInfoUsed$start)

        colorsAll = rep(CPG_SEA_COLOR, NROW(validList))
        colorsAll[meth[validList,"start"] <= cpgInfoUsed$end[cpgId] + 4000 | meth[validList,"end"] >= cpgInfoUsed$start[cpgId+1] - 4000 ] =
                CPG_SHELVE_COLOR
        colorsAll[meth[validList,"start"] <= cpgInfoUsed$end[cpgId] + 2000 | meth[validList,"end"] >= cpgInfoUsed$start[cpgId+1] - 2000 ] =
                CPG_SHORE_COLOR
        colorsAll[cpgInfoUsed$end[cpgId] >= meth[validList,"start"]] = CPG_ISLAND_COLOR
    } else {
        colorsAll = rep("red", NROW(validList))
    }

    samp = sample(1:NROW(groupAll),NROW(groupAll))
    lineType = ifelse(noDataLine == TRUE, "p", "b")

    if(!tooWide)
    for(j in 1:NROW(groupAll)) {
        i = samp[j]

        coverage   = meth@data.numCs[validList,groupAll[i]] + meth@data.numTs[validList,groupAll[i]]
        cNums      = meth@data.numCs[validList,groupAll[i]]
        start      = meth@data.start[validList]
        nonNAdata = which(coverage>0)
        coverage = coverage[nonNAdata]
        cNums = cNums[nonNAdata]
        start = start[nonNAdata]
        colors = colorsAll[nonNAdata]
        ord = order(start)
        if(i <= NROW(group1)) {
            points(start[ord], (cNums/coverage*100)[ord], type=lineType, pch=pchList[i], col=colors[ord], cex=sqrt(coverage)/10*cex, lty=2)
        } else points(start[ord], (cNums/coverage*100)[ord], type=lineType, pch=pchList[i], col="black", cex=sqrt(coverage)/10*cex, lty=2)
    }

    if(!missing(cpgInfo)) {
        CpG_shelf = 4000
        CpG_shore = 2000
        cpgInfoToDraw = which(cpgInfo[[2]]$chr == chr & cpgInfo[[2]]$end + CpG_shelf > loc.range[1] & cpgInfo[[2]]$start - CpG_shelf < loc.range[2])
        rect(xleft=loc.range[1],xright=loc.range[2], ybottom=yLLoc, ytop=yLLoc*0.75+yHLoc*0.25, col="blue", border="blue")
        if(length(cpgInfoToDraw) > 0) {
            rect(xleft=c(cpgInfo[[2]]$start[cpgInfoToDraw]-CpG_shelf,cpgInfo[[2]]$end[cpgInfoToDraw]+CpG_shore),
                      xright=c(cpgInfo[[2]]$start[cpgInfoToDraw]-CpG_shore,cpgInfo[[2]]$end[cpgInfoToDraw]+CpG_shelf),
                      ybottom=yLLoc, ytop=(yLLoc+yHLoc)*0.5, col=CPG_SHELVE_COLOR, border=CPG_SHELVE_COLOR)
            rect(xleft=c(cpgInfo[[2]]$start[cpgInfoToDraw]-CpG_shore, cpgInfo[[2]]$end[cpgInfoToDraw]),
                      xright=c(cpgInfo[[2]]$start[cpgInfoToDraw], cpgInfo[[2]]$end[cpgInfoToDraw]+CpG_shore),
                      ybottom=yLLoc, ytop=yLLoc*0.25+yHLoc*0.75, col=CPG_SHORE_COLOR, border=CPG_SHORE_COLOR)
            rect(xleft=cpgInfo[[2]]$start[cpgInfoToDraw], xright=cpgInfo[[2]]$end[cpgInfoToDraw],
                      ybottom=yLLoc, ytop=yHLoc, col=CPG_ISLAND_COLOR, border=CPG_ISLAND_COLOR)
        }

        axis(2, at=(yHLoc+yLLoc)/2, labels = "CpG", las=1, tick=FALSE)
        yHLoc = yHLoc-EXTRA
        yLLoc = yLLoc-EXTRA
    }

    if(!missing(refGeneInfo)) {
        whichChr = which(refGeneInfo[[1]]$chr == chr)
        whichRange = refGeneInfo[[1]]$start[whichChr]:refGeneInfo[[1]]$end[whichChr]
        promoter_length = 2000
        COLOR_INTERGENIC = "gray47"
        ### intergenic
        refseqToDraw = whichRange[refGeneInfo[[2]]$end[whichRange] > (loc.range[1] - promoter_length) & refGeneInfo[[2]]$start[whichRange] < (loc.range[2] + promoter_length)]
        rect(xleft=loc.range[1],xright=loc.range[2], ybottom=yLLoc, ytop=yLLoc+(yHLoc-yLLoc)/4,
                  col=COLOR_INTERGENIC, border=COLOR_INTERGENIC)

        if(length(refseqToDraw) > 0) {
            whichStrandPlus = (refGeneInfo[[2]]$strand[refseqToDraw] == "+")

            ###exon
            exon_start = refGeneInfo[[2]]$exon.index.start[refseqToDraw[1]]
            exon_end   = refGeneInfo[[2]]$exon.index.start[refseqToDraw[NROW(refseqToDraw)]] + refGeneInfo[[2]]$exon.num[refseqToDraw[NROW(refseqToDraw)]] - 1
            rect(xleft=refGeneInfo[[3]]$exon_start[exon_start:exon_end],
                 xright=refGeneInfo[[3]]$exon_end[exon_start:exon_end],
                 ybottom=yLLoc, ytop=yHLoc, col="brown1", border="brown1")

            ### 3' UTR
            utrList = c(!whichStrandPlus&(refGeneInfo[[2]]$exon.num[refseqToDraw]>0), whichStrandPlus&(refGeneInfo[[2]]$exon.num[refseqToDraw]>0))

            if(sum(utrList) > 0)
            rect(xleft=c(refGeneInfo[[2]]$start[refseqToDraw], refGeneInfo[[2]]$cds_end[refseqToDraw])[utrList],
                xright=c(refGeneInfo[[2]]$cds_start[refseqToDraw], refGeneInfo[[2]]$end[refseqToDraw])[utrList],
                ybottom=yLLoc, ytop=(yLLoc+yHLoc)*0.5, col="yellow3", border="yellow3")

            ### promoter
            rect(xleft=c(refGeneInfo[[2]]$start[refseqToDraw[whichStrandPlus]]-promoter_length, refGeneInfo[[2]]$end[refseqToDraw[!whichStrandPlus]]),
                     xright=c(refGeneInfo[[2]]$start[refseqToDraw[whichStrandPlus]], refGeneInfo[[2]]$end[refseqToDraw[!whichStrandPlus]]+promoter_length),
                     ybottom=yLLoc, ytop=(yLLoc+yHLoc*3)/4, col="green", border="green")

            ### noncoding RNA
            utrList = refGeneInfo[[2]]$exon.num[refseqToDraw]==0

            if(sum(utrList) > 0)
            rect(xleft=refGeneInfo[[2]]$start[refseqToDraw][utrList],
                xright=refGeneInfo[[2]]$end[refseqToDraw][utrList],
                ybottom=yLLoc, ytop=(yLLoc+yHLoc)*0.5, col="yellow4", border="yellow4")


            ### 5' UTR
            utrList = c(whichStrandPlus&(refGeneInfo[[2]]$exon.num[refseqToDraw]>0), !whichStrandPlus&(refGeneInfo[[2]]$exon.num[refseqToDraw]>0))

            if(sum(utrList) > 0)
            rect(xleft=c(refGeneInfo[[2]]$start[refseqToDraw], refGeneInfo[[2]]$cds_end[refseqToDraw])[utrList],
                xright=c(refGeneInfo[[2]]$cds_start[refseqToDraw], refGeneInfo[[2]]$end[refseqToDraw])[utrList],
                ybottom=yLLoc, ytop=(yLLoc+yHLoc)*0.5, col="yellow", border="yellow")

             ###CDS
            rect(xleft=refGeneInfo[[2]]$cds_start[refseqToDraw], xright=refGeneInfo[[2]]$cds_end[refseqToDraw], ybottom=yLLoc, ytop=(yLLoc+yHLoc)*0.5, col="lightskyblue", border="lightskyblue")
        }
        axis(2, at=(yHLoc+yLLoc)/2, labels = "refGene", las=1, tick=FALSE)

        yHLoc = yHLoc-EXTRA
        yLLoc = yLLoc-EXTRA
    }

    if(!missing(myDiff)) {
        abline(h=c(yHLoc, yLLoc))
        validList <- which((myDiff[,"chr"] == chr) & (myDiff[,"start"] >= loc.range[1]) & (myDiff[,"end"] <= loc.range[2]))
        if(length(validList)>0) {
        #### Draw oly when there are data here in this range
            validList = validList[order(myDiff[validList,"start"])]
            if(noGroupEst == FALSE) {
                muName = paste("mu", groups, sep="")
                lines(myDiff[validList,"start"], myDiff[validList,muName[1]], type="b", pch=16, col=rgb(255,0,0,100,maxColorValue=255), cex=cex, lwd=1, lty=1)
                lines(myDiff[validList,"start"], myDiff[validList,muName[2]], type="b", pch=16, col=rgb(0,0,0,100,maxColorValue=255), cex=cex, lwd=1, lty=1)
            }
            colorAll = rep("blue", NROW(validList))
            colorAll[myDiff[validList, "qvalue"]<sigQ] = "red"
            yBottom = rep(yLLoc+1.25, NROW(validList))
            yTop = yLLoc - log10(pmax(myDiff[validList, "pvalue"], 1e-10))/10*(yHLoc-yLLoc)
            rect(xleft=myDiff[validList,"start"], xright=myDiff[validList,"start"], ybottom=yLLoc, ytop=yTop, col=colorAll, border=colorAll)
        }
        axis(2, at=(yHLoc+yLLoc)/2, labels = expression(-log[10](pvalue)), las=1, tick=FALSE)
        mtext(text=c("10","0"), side=2, at=c(yHLoc,yLLoc), line=0.1, cex.axis=0.8, las=1)

        yHLoc = yHLoc-EXTRA
        yLLoc = yLLoc-EXTRA
    }

    if(!missing(tfbsInfo)) {
        whichChr = which(tfbsInfo[[1]]$chr == chr)
        whichRange = tfbsInfo[[1]]$start[whichChr]:tfbsInfo[[1]]$end[whichChr]

        tfbsId = whichRange[loc.range[1] <= tfbsInfo[[2]]$chromStart[whichRange] &loc.range[2] >= tfbsInfo[[2]]$chromStart[whichRange]]

        if(length(tfbsId) > 0) {
             if(tfbsDense==FALSE) {
                 colorStrand = rep(rgb(0,0,0,100,"",255), NROW(tfbsId))
                 colorStrand[tfbsInfo[[2]]$strand[tfbsId] == "-"] = rgb(0,0,255,100,"",255)
                 colorStrand[tfbsInfo[[2]]$strand[tfbsId] == "+"] = rgb(255,155,0,100,"",255)
                 rect(xleft=tfbsInfo[[2]]$chromStart[tfbsId], xright=tfbsInfo[[2]]$chromEnd[tfbsId],
                      ybottom=yHLoc-(0:(NROW(tfbsId)-1))%%10 * 3 - 3, ytop=yHLoc-(0:(NROW(tfbsId)-1))%%10 * 3-1,
                      col=colorStrand, border="transparent")
                 text(tfbsInfo[[2]]$chromStart[tfbsId], yHLoc-(0:(NROW(tfbsId)-1))%%10 * 3 - 2, tfbsInfo[[2]]$name[tfbsId], cex=0.6, pos=2, offset=0)
              } else {
                 colorStrand = rep(rgb(0,0,0,40,"",255), NROW(tfbsId))
                 colorStrand[tfbsInfo[[2]]$strand[tfbsId] == "-"] = rgb(0,0,255,40,"",255)
                 colorStrand[tfbsInfo[[2]]$strand[tfbsId] == "+"] = rgb(255,155,0,40,"",255)

                 rect(xleft=tfbsInfo[[2]]$chromStart[tfbsId], xright=tfbsInfo[[2]]$chromEnd[tfbsId],
                      ybottom=yLLoc, ytop=yHLoc, col=colorStrand, border=colorStrand)
               }
#                  col=colorStrand[tfbsInfo[[2]]$strand[tfbsId]], border="red")
        }

        axis(2, at=(yHLoc+yLLoc)/2, labels = "tfbs", las=1, tick=FALSE)
    }

    title(xlab = paste("base(",chr,")", sep=""))
}

#' Read RefGene model information.
#'
#' This function reads RefGene model information for future annotation.
#'
#'   This function reads RefGene informaiton for future annotation. Same genome assembly of data should be used in annotation. This function does not check genome assembly type nor version.
#'
#' @param infile A path to file that contains refGene information.
#'
#' @return RefGene information object that can be used in function \code{\link{refGeneAnnotation}}.
#'
#' @seealso \code{\link{refGeneAnnotation}} and \code{\link{refGeneAnnotationPlot}}
#'
#' @examples
#' refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt",
#'                             package = "methylSig"))
#'
#' @keywords refGeneAnnotation
#'
#' @export
getRefgeneInfo <- function(infile) {
    refgen <- read.table(file=infile, sep="\t")
    refgen = refgen[order(refgen$V5+as.numeric(refgen$V3)*300000000),]

## Noncoding exon is useful
##
#    nonCodingList = which(refgen$V7 == refgen$V8)
#    refgen$V9[nonCodingList] = 0
#    refgen$V10[nonCodingList] = refgen$V11[nonCodingList] = NA

    refGeneInfo = list()
    refGeneInfo[[1]] = data.frame(chr = unique(refgen$V3))
    refGeneInfo[[1]]$start = match(refGeneInfo[[1]]$chr, refgen$V3)
    refGeneInfo[[1]]$end = c(refGeneInfo[[1]]$start[-1]-1,NROW(refgen))

    refGeneInfo[[2]] = data.frame(chr=refgen$V3, strand=refgen$V4, start=refgen$V5, end=refgen$V6, cds_start=refgen$V7,
                                   cds_end=refgen$V8, exon.num=refgen$V9, name=as.character(refgen$V2), gene=as.character(refgen$V13))
    refGeneInfo[[3]] = data.frame(exon_start = as.numeric(do.call(c, strsplit(as.character(refgen$V10[refgen$V9!=0]),","))),
                                   exon_end = as.numeric(do.call(c, strsplit(as.character(refgen$V11[refgen$V9!=0]),","))))

    refGeneInfo[[2]]$exon.index.start = cumsum(c(1,refGeneInfo[[2]]$exon.num[1:(NROW(refGeneInfo[[2]]$exon.num)-1)]))

    refGeneInfo
}

#' Read transcription factor binding site information
#'
#' This function reads the transcription factor binding site information for future annotation.
#'
#' @param infile A path to BED4+2 file (\url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}) containing transcription factor binding site information.
#'
#' @return An object containing transcription factor binding sites information is returned. It can be used in function \code{\link{methylSig.tfbsEnrichTest}}.
#'
#' @seealso \code{\link{methylSig.tfbsEnrichTest}}
#'
#' @examples
#'    tfbsInfo = getTFBSInfo(system.file("annotation", "tfbsUniform.txt",
#'                        package = "methylSig"))
#'
#' @keywords tfbsAnnotation
#'
#' @export
getTFBSInfo <- function(infile) {
    tfbs <- read.table(file=infile)

    if(is.numeric(tfbs[1,1])) {
        tfbs = tfbs[,2:7,]
    } else {
        tfbs = tfbs[,1:6]
    }
#    names(tfbs) = c("bin", "chrom", "chromStart", "chromEnd", "name", "score", "strand", "zScore")

    names(tfbs) = c("chrom", "chromStart", "chromEnd", "name", "score", "strand")
#    tfbs$name=sub("^..","",tfbs$name)

    tfbsInfo = list()
    tfbsInfo[[1]] = data.frame(chr = unique(tfbs$chrom))
    tfbsInfo[[1]]$start = match(tfbsInfo [[1]]$chr, tfbs$chrom)
    tfbsInfo[[1]]$end = c(tfbsInfo[[1]]$start[-1]-1,NROW(tfbs))

    tfbsInfo[[2]] = tfbs
    if(class(tfbsInfo[[2]]$name) == "character") tfbsInfo[[2]]$name=as.factor(tfbsInfo[[2]]$name)
    tfbsInfo
}

#' Read CpG island information from a text file.
#'
#' This function reads CpG island information from a text or bed file for future annotation.
#'
#' @param infile A path to a BED file containing chromosome, start, and end positions for CpG islands.
#'
#' @return CpG island information database that can be used in function cpgAnnotation.
#'
#' @section Warning: The same genome assembly of data should be used in annotation. This function does not check genome assembly type or version.
#'
#' @seealso \code{\link{cpgAnnotationPlot}} and \code{\link{cpgAnnotation}}
#'
#' @examples
#' cpgInfo = getCpGInfo(system.file("annotation", "cpgi.hg18.bed.txt",
#'                       package = "methylSig"))
#'
#' @keywords CpGAnnotation
#'
#' @export
getCpGInfo <-function(infile) {
    CpGInfo= list()
    input = read.table(file=infile)
    startCol = ifelse(is.numeric(input[1,1]), 2,1)
    CpGInfo[[2]] <- data.frame(chr=input[,startCol], start=input[,startCol+1], end=input[,startCol+2])

    CpGInfo[[1]] = data.frame(chr = unique(CpGInfo[[2]]$chr))
    CpGInfo[[1]]$start = match(CpGInfo[[1]]$chr, CpGInfo[[2]]$chr)
    ####### CpG island exact ending site
    CpGInfo[[1]]$end = c(CpGInfo[[1]]$start[-1]-1,NROW(CpGInfo[[2]]))-1

    CpGInfo
}

# Called by methylSig.tfbsEnrichTest
getTFBSCountByChrom <- function(tfbsInfo, startEnd, listToCount) {
    maxOverLaps = 200
    cpgId = findInterval(listToCount,tfbsInfo[[2]]$chromStart[startEnd])
    cpgIdAll = pmax(0,rep(cpgId,each=maxOverLaps)+rep(-(maxOverLaps-1):0,NROW(cpgId)))
    whichValid = which(rep(listToCount,each=maxOverLaps)<=c(0,tfbsInfo[[2]]$chromEnd[startEnd])[cpgIdAll+1])

    #### In case of that not all slots have data
    table(c(1:NROW(levels(tfbsInfo[[2]]$name)),tfbsInfo[[2]]$name[startEnd][cpgIdAll[whichValid]]))-1
}

# Called by methylSig.tfbsEnrichTest
is.TFBSByChrom <- function(tfbsInfo, startEnd, listToCount) {
    maxOverLaps = 200
    cpgId = findInterval(listToCount,tfbsInfo[[2]]$chromStart[startEnd])
    cpgIdAll = pmax(0,rep(cpgId,each=maxOverLaps)+rep(-(maxOverLaps-1):0,NROW(cpgId)))
    whichValid = rep(listToCount,each=maxOverLaps)<=c(0,tfbsInfo[[2]]$chromEnd[startEnd])[cpgIdAll+1]

    #### In case of that not all slots have data
    colSums(matrix(whichValid, nrow=maxOverLaps))>0
}

#' CpG island Annotation pie chart
#'
#' This function generates a pie chart that shows the proportions of CpG sites annotated to a CpG island, CpG shore, CpG shelf or interCGI region.
#'
#' CpG shore is defined as the region outside of a CpG island but within 2,000 bp of any CpG island. CpG shelf is defined as the region 2,000-4,000bp away from a CPG island. When overlaped, the priority is CpG island, then CpG shore.
#'
#' @param listFrom An R object generated by \code{\link{cpgAnnotation}} function.
#' @param main An overall title for the plot. The default is "Plot".
#' @param color Colors used in the plot to show CpG island, CpG shore, CpG shelve and interCGI regions. If noShelf = TRUE, then only three colors are needed.
#' @param noShelf Whether CpG shelves should be included in the plot. If noShelf = TRUE, then CpG shelves are treated as part of interCGI region. Default is FALSE.
#'
#' @seealso \code{\link{cpgAnnotation}} and \code{\link{getCpGInfo}}
#'
#' @examples
#' data(sampleData)
#'
#' cpgInfo = getCpGInfo(system.file("annotation", "cpgi.hg18.bed.txt",
#'                        package = "methylSig"))
#'
#' cpgAnn = cpgAnnotation(cpgInfo,myDiffSigboth)
#'
#' cpgAnnotationPlot(cpgAnn,main="ALL")
#'
#' @keywords CpGAnnotation
#'
#' @export
cpgAnnotationPlot <- function(listFrom, main="Plot", color, noShelf = FALSE) {
   CPG_ISLAND_COLOR="chartreuse4"
   CPG_SHORE_COLOR="green3"
   CPG_SHELVE_COLOR="gold2"
   CPG_SEA_COLOR="BLUE"

   if(noShelf) {
        slices <- c(sum(listFrom== "CpGIsland"),
                    sum(listFrom== "CpGShore"),
                    sum(listFrom== "CpGShelf") + sum(listFrom== "InterCGI"))
        lbls <- c("CpG Island", "CpG Shore", "InterCGI")
        if(missing(color)) color = c(CPG_ISLAND_COLOR,CPG_SHORE_COLOR,CPG_SEA_COLOR="BLUE")
    } else {
        slices <- c(sum(listFrom== "CpGIsland"),
                    sum(listFrom== "CpGShore"),
                    sum(listFrom== "CpGShelf"),
                    sum(listFrom== "InterCGI"))
        lbls <- c("CpG Island", "CpG Shore", "CpG Shelf", "InterCGI")
        if(missing(color)) color = c(CPG_ISLAND_COLOR,CPG_SHORE_COLOR,CPG_SHELVE_COLOR,CPG_SEA_COLOR="BLUE")
    }
    pct <- round(slices/sum(slices)*100)
    lbls <- paste(lbls, pct) # add percents to labels
    lbls <- paste(lbls,"%",sep="") # ad % to labels
    pie(slices,labels = lbls, col=color, main=main)
}

#' CpG island annotation function
#'
#' This function annotates the CpG sites or regions to CpG island, CpG shore, CpG shelf or interCGI regions.
#'
#'   This function generates an object to annotate CpG island information. The same genome assembly should be used for `myDiff' and `cpgInfo'. CpG shores are defined as the region outside CpG islands but within 2,000 bp of any CpG island. CpG shelves are defined as the region within 2,000 bp away from a CpG shore. When regions overlap, the priority is CpG island, then CpG shore. Regions >4,000 bp from a CpG island are deemed interCGI.
#'
#' @param cpgInfo A CpG island information object.
#' @param myDiff A `methylSigDiff' object to be used for annotation.
#'
#' @return A factor vector with the same size as `myDiff' that indicates whether the related sites in `myDiff' are in a CpG island, CpG shore, CpG shelf or interCGI region.
#'
#' @seealso \code{\link{cpgAnnotationPlot}} and \code{\link{getCpGInfo}}
#'
#' @examples
#' data(sampleData)
#'
#' cpgInfo = getCpGInfo(system.file("annotation", "cpgi.hg18.bed.txt",
#'                                  package = "methylSig"))
#'
#' cpgAnn = cpgAnnotation(cpgInfo,myDiffSigboth)
#'
#' @keywords CpGAnnotation
#'
#' @export
cpgAnnotation <- function(cpgInfo, myDiff) {
    retList = rep(factor(c("Unknown", "CpGIsland","CpGShore", "CpGShelf", "InterCGI"))[1],NROW(myDiff@data.ids))
    for(chr in levels(myDiff@data.chr)) {
        whichChr = which(chr == cpgInfo[[1]]$chr)
        if(length(whichChr) > 0) {
            listChr = which(myDiff@data.chr == chr)
            retList[listChr] = cpgAnnotationByChrom(cpgInfo, cpgInfo[[1]]$start[whichChr]:cpgInfo[[1]]$end[whichChr], myDiff@data.start[listChr])
        }
    }
    retList
}

# Called by cpgAnnotation
cpgAnnotationByChrom <- function(cpgInfo, startEnd, listToAnnotat) {
    CpG_shelf = 4000
    CpG_shore = 2000

    ret <- rep(factor(c("Unknown", "CpGIsland","CpGShore", "CpGShelf", "InterCGI"))[5], NROW(listToAnnotat))
    cpgInfoUsed = data.frame(start = c(-Inf, cpgInfo[[2]]$start[startEnd], Inf), end = c(-Inf, cpgInfo[[2]]$end[startEnd], Inf))
    cpgId = findInterval(listToAnnotat, cpgInfoUsed$start)

    ret[listToAnnotat <= cpgInfoUsed$end[cpgId] + 4000 | listToAnnotat >= cpgInfoUsed$start[cpgId+1] - 4000 ] = "CpGShelf"
    ret[listToAnnotat <= cpgInfoUsed$end[cpgId] + 2000 | listToAnnotat >= cpgInfoUsed$start[cpgId+1] - 2000 ] = "CpGShore"
    ret[cpgInfoUsed$end[cpgId] >= listToAnnotat] = "CpGIsland"
    ret
}

#' RefGene annotation plot
#'
#' This function generates a refGene annotation bar chart.
#'
#' Promoter is defined as the region 1,000bp upstream from a transcription start site (TSS) to the TSS. When regions overlap, the priority to assign regions for the CpG sites is defined by `priority' argument.
#'
#' @param listFrom Object from fuction `refGeneAnnotation'.
#' @param main An overall title for the plot. Default is "Plot".
#' @param priority A character vector indicating the priority order for assigning the annotated CpG site when overlapping occurs.  Default is coding sequence exon, promoter, noncoding exon, 5' UTR exon, 3' UTR exon, intron and then integenic region.  Only "cds", "promoter", "noncoding", "5'utr" and "3'utr" are acceptable.
#'
#' @seealso \code{\link{refGeneAnnotation}} and \code{\link{getRefgeneInfo}}
#'
#' @examples
#' data(sampleData)
#'
#' refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt",
#'                  package = "methylSig"))
#'
#' refGeneAnn = refGeneAnnotation(refGeneInfo, myDiffSigboth)
#'
#' refGeneAnnotationPlot(refGeneAnn,main="ALL", priority=c("promoter","cds",
#'                            "noncoding", "5'utr", "3'utr"))
#'
#' @keywords refGeneAnnotation
#'
#' @export
refGeneAnnotationPlot <- function(listFrom, main="Plot", priority=c("cds", "promoter","noncoding", "5'utr", "3'utr")) {
    countInfo= rep(NA, NCOL(listFrom))

    INTERGENIC = 0
    INTRON     = 1
    UTR3       = 2
    UTR5       = 3
    NONCODING  = 4
    CDS        = 5
    PROMOTER   = 6

    countInfo[colSums(listFrom)==0] = INTERGENIC
    countInfo[colSums(listFrom)> 0 & (listFrom["exon",]==FALSE)] = INTRON

    for(whatNext in 5:1) {
        if(priority[whatNext] == "promoter") {
           countInfo[which(listFrom["promoter",])] = PROMOTER
        } else if(priority[whatNext] == "3'utr") {
            countInfo[listFrom["utr3",] & listFrom["exon",]] = UTR3
        } else if(priority[whatNext] == "cds") {
            countInfo[listFrom["cds",] & listFrom["exon",]] = CDS
        } else if(priority[whatNext] == "5'utr") {
            countInfo[listFrom["utr5",] & listFrom["exon",]] = UTR5
        } else if(priority[whatNext] == "noncoding") {
            countInfo[listFrom["noncoding",] & listFrom["exon",]] = NONCODING
        } else {
            stop("priority includes \"cds\", \"promoter\", \"noncoding\", \"utr5\" and \"utr3\" only")
        }
    }

    slices <- c(sum(countInfo==INTERGENIC,na.rm=TRUE),sum(countInfo==INTRON, na.rm=TRUE)) ### integenic and intron
    lbls <- c("Intergenic", "Intron")
    for(whatNext in 5:1) {
        if(priority[whatNext] == "promoter") {
            slices = c(slices, sum(countInfo==PROMOTER,na.rm=TRUE))
            lbls = c(lbls, "Promoter")
        } else if(priority[whatNext] == "3'utr") {
            slices = c(slices, sum(countInfo==UTR3,na.rm=TRUE))
            lbls = c(lbls, "3'UTR")
        } else if(priority[whatNext] == "cds") {
            slices = c(slices, sum(countInfo==CDS,na.rm=TRUE))
            lbls = c(lbls, "CDS")
        } else if(priority[whatNext] == "5'utr") {
            slices = c(slices, sum(countInfo==UTR5,na.rm=TRUE))
            lbls = c(lbls, "5'UTR")
        } else if(priority[whatNext] == "noncoding") {
            slices = c(slices, sum(countInfo==NONCODING,na.rm=TRUE))
            lbls = c(lbls, "Nonconding")
        }
    }

    pct <- paste(round(slices/sum(slices)*100),"%",sep="")
    par(las=2)
    b <- barplot(slices,col=rainbow(length(lbls)), main=main, ylim=c(0,max(slices)*1.1), xaxt = "n")
    text(x=b,y=slices,labels=pct,pos=3)
    text(x=b, y=-max(slices)*0.02, lbls, xpd=TRUE, srt=45, pos=2)
}



#' Annotation function using RefGene models
#'
#' This function annotates the CpG sites to promoter, coding DNA sequence (CDS), exons, 3' untranslated region, 5' untranslated region, nocoding RNA region, or integenic regions using the RefGene gene models.
#'
#'   This function generates an object to annotate refGene information. The same genome assembly should be used for `myDiff' and `refGeneInfo'. The promoter region is defined as 1,000bp upstream from a transcription starting site (TSS) to the TSS.
#'
#' @param refGeneInfo RefGene information object from function \code{\link{getRefgeneInfo}}.
#' @param myDiff `methylSigDiff' object to be used for annotation.
#'
#' @return A logical matrix of 6 rows and number of columns equal to the size of `myDiff' with each entry indicating whether the related sites in `myDiff' are in a promoter, CDS, 3' UTR, 5' UTR, noncoding RNA or exon regions. Regions can overlap because of a single gene can have multiple isoforms.
#'
#' @seealso \code{\link{refGeneAnnotationPlot}} and \code{\link{getRefgeneInfo}}
#'
#' @examples
#' data(sampleData)
#' refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt",
#'                                package = "methylSig"))
#'
#' refGeneAnn = refGeneAnnotation(refGeneInfo,myDiffSigboth)
#'
#' @keywords refGeneAnnotation
#'
#' @export
refGeneAnnotation <- function(refGeneInfo, myDiff) {
    refInfoCountInfo = matrix(NA,nrow=6, ncol=NROW(myDiff@data.chr))

    for(chr in levels(myDiff@data.chr)) {
        vlist = which(myDiff@data.chr==chr)

        refChrId = which(refGeneInfo[[1]]$chr == chr)
        if(length(refChrId) > 0) {
            startEnd = refGeneInfo[[1]]$start[refChrId]: refGeneInfo[[1]]$end[refChrId]
            refInfoCountInfo[,vlist] = refGeneAnnotationByChrom(refGeneInfo, startEnd, myDiff@data.start[vlist])
        }
    }

    rownames(refInfoCountInfo) = c("promoter", "cds", "utr3", "utr5", "exon", "noncoding")
    refInfoCountInfo
}

# Called by refGeneAnnotation
refGeneAnnotationByChrom <- function(refGeneInfo, startEnd, listToAnnotat) {
    ### intergenic
    ret <- matrix(FALSE, nrow=6, ncol=NROW(listToAnnotat))
    rownames(ret) = c("promoter", "cds", "utr3", "utr5", "exon", "noncoding")
    ##
    ordList = order(listToAnnotat)
    listToAnnotat= listToAnnotat[ordList]
    promoter_length = 2000

    startMatch = findInterval(refGeneInfo[[2]]$start[startEnd], listToAnnotat)
    endMatch = findInterval(refGeneInfo[[2]]$end[startEnd], listToAnnotat)
    cdsstartMatch = findInterval(refGeneInfo[[2]]$cds_start[startEnd], listToAnnotat)
    cdsendMatch = findInterval(refGeneInfo[[2]]$cds_end[startEnd], listToAnnotat)
    promoterStart = refGeneInfo[[2]]$cds_start[startEnd] - promoter_length
    reverseStrand = refGeneInfo[[2]]$strand[startEnd] == "-"
    promoterStart[reverseStrand] = refGeneInfo[[2]]$end[startEnd][reverseStrand] + promoter_length

    promoterMatch = findInterval(promoterStart, listToAnnotat)

    #### Promoter

    promoterList =  !reverseStrand & (startMatch > promoterMatch)
    for(i in which(promoterList)) {
        ret["promoter",(promoterMatch[i]+1):startMatch[i]] = TRUE
    }
    promoterList =  reverseStrand & (endMatch < promoterMatch)
    for(i in which(promoterList)) {
        ret["promoter",(endMatch[i]+1):promoterMatch[i]] = TRUE
    }

    for(i in which(cdsendMatch > cdsstartMatch)) {
        ret["cds",(cdsstartMatch[i]+1):cdsendMatch[i]] = TRUE
    }

    #### Noncoding

    for(i in which(endMatch > startMatch & (refGeneInfo[[2]]$cds_start[startEnd]==refGeneInfo[[2]]$cds_end[startEnd]))) {
        ret["noncoding",(startMatch[i]+1):endMatch[i]] = TRUE
    }

    ### 5' UTR
    utrList =  !reverseStrand & (cdsstartMatch > startMatch) & (refGeneInfo[[2]]$cds_start[startEnd]!=refGeneInfo[[2]]$cds_end[startEnd])
    for(i in which(utrList)) {
        ret["utr5",(startMatch[i]+1):cdsstartMatch [i]] = TRUE
    }
    utrList =  reverseStrand & (endMatch > cdsendMatch) & (refGeneInfo[[2]]$cds_start[startEnd]!=refGeneInfo[[2]]$cds_end[startEnd])
    for(i in which(utrList)) {
        ret["utr5",(cdsendMatch[i]+1):endMatch [i]] = TRUE
    }

    ### 3' UTR
    utrList =  reverseStrand & (cdsstartMatch > startMatch) & (refGeneInfo[[2]]$cds_start[startEnd]!=refGeneInfo[[2]]$cds_end[startEnd])
    for(i in which(utrList)) {
        ret["utr3",(startMatch[i]+1):cdsstartMatch [i]] = TRUE
    }
    utrList =  !reverseStrand & (endMatch > cdsendMatch) & (refGeneInfo[[2]]$cds_start[startEnd]!=refGeneInfo[[2]]$cds_end[startEnd])
    for(i in which(utrList)) {
        ret["utr3",(cdsendMatch[i]+1):endMatch [i]] = TRUE
    }

    #### Exon #####
    lastOfStartEnd = startEnd[NROW(startEnd)]
    validListInExon = refGeneInfo[[2]]$exon.index.start[startEnd[1]]:( refGeneInfo[[2]]$exon.index.start[lastOfStartEnd]
                                                                     +refGeneInfo[[2]]$exon.num        [lastOfStartEnd]-1)
    startMatch = findInterval(refGeneInfo[[3]]$exon_start[validListInExon], listToAnnotat)
    endMatch   = findInterval(refGeneInfo[[3]]$exon_end[validListInExon], listToAnnotat)

    exonList =  (endMatch > startMatch)
    for(i in which(exonList)) {
        ret["exon",(startMatch[i]+1):endMatch[i]] = TRUE
    }

    ret[,ordList] = ret
    ret
}

#' Perform transcription factor enrichment test among differentially methylated cytosines or regions
#'
#' This function tests for enriched transcription binding sites among differentially methylated sites or regions using a binomial test. It also has an option to generate a bar plot to show the significance levels for the enriched transcription factor binding sites.
#'
#' Likelihood ratio test is used based on the binomial distribution.
#'
#' @param myDiff methylSigDiff object that contains all CpG sites that are tested for differential methylation.
#' @param dmcList A boolean vector of the same length as myDiff defining the DMCs or DMRs. A value of TRUE for dmcList[i] means that the ith CpG site in myDiff is a DMC.
#' @param tfbsInfo An object that contains transcription factor binding sites information.
#' @param plot If plot = TRUE, this function draws a plot. Otherwise no plot is given. Default is TRUE.
#' @param max.plot.num Maximum number of transcription factor binding sites on the plot. Default is 10.
#' @param main An overall title for the plot.
#'
#' @return p-values of enriched transcription factor binding sites.
#'
#' @seealso \code{\link{getTFBSInfo}}
#'
#' @examples
#' data(sampleData)
#'
#' tfbsInfo = getTFBSInfo(system.file("annotation", "tfbsUniform.txt",
#'                         package = "methylSig"))
#'
#' pvalue = methylSig.tfbsEnrichTest(myDiffSigboth,
#'          myDiffSigboth[,"qvalue"]<0.05 & abs(myDiffSigboth[,"meth.diff"])>25,
#'          tfbsInfo)
#'
#' @keywords tfbsAnnotation
#'
#' @export
methylSig.tfbsEnrichTest <- function(myDiff, dmcList, tfbsInfo, plot=TRUE, max.plot.num=10, main="tfbs Enrichment") {
    tfbsInfoCountInfo = NULL
    tfbsInfoCountInfoDMC = NULL
    tfbsInfoIsTfbsDMC = NULL
    tfbsInfoIsTfbsTotal = NULL

    for(chr in levels(myDiff@data.chr)) {
        tfbsChrId = which(tfbsInfo[[1]]$chr == as.character(chr))
        if(length(tfbsChrId) > 0) {
            vlist = (myDiff@data.chr==chr)
            startEnd = tfbsInfo[[1]]$start[tfbsChrId]: tfbsInfo[[1]]$end[tfbsChrId]
            tfbsInfoCountInfo = cbind(tfbsInfoCountInfo,getTFBSCountByChrom(tfbsInfo, startEnd, myDiff@data.start[vlist]))

            tfbsInfoIsTfbsTotal = c(tfbsInfoIsTfbsTotal, is.TFBSByChrom(tfbsInfo, startEnd, myDiff@data.start[vlist]))

            vlist = (myDiff@data.chr==chr & dmcList)
            startEnd = tfbsInfo[[1]]$start[tfbsChrId]: tfbsInfo[[1]]$end[tfbsChrId]
            tfbsInfoCountInfoDMC= cbind(tfbsInfoCountInfoDMC,getTFBSCountByChrom(tfbsInfo, startEnd, myDiff@data.start[vlist]))

            tfbsInfoIsTfbsDMC = c(tfbsInfoIsTfbsDMC, is.TFBSByChrom(tfbsInfo, startEnd, myDiff@data.start[vlist]))
        }
    }

    tfbsInfoCountInfoSum = apply(tfbsInfoCountInfo, 1, sum)
    tfbsInfoCountInfoDMCSum = apply(tfbsInfoCountInfoDMC, 1, sum)
    totalALL = sum(tfbsInfoIsTfbsTotal)
    totalDMC =  sum(tfbsInfoIsTfbsDMC)

    rates = tfbsInfoCountInfoSum/totalALL
    logLik = 2*(tfbsInfoCountInfoDMCSum*
                  log(pmax(tfbsInfoCountInfoDMCSum,1e-100)/totalDMC/rates)
          + (totalDMC - tfbsInfoCountInfoDMCSum)*log((1-tfbsInfoCountInfoDMCSum/totalDMC)/(1-rates)))

    pvalue = pchisq(logLik, 1, lower.tail=FALSE)
    pvalue [tfbsInfoCountInfoDMCSum/totalDMC<tfbsInfoCountInfoSum/totalALL] = 1

    names(pvalue) = levels(tfbsInfo[[2]]$name)
    pvalue = pvalue[tfbsInfoCountInfoSum > 0]

    numPlot = min(10, sum(pvalue < 1))
    pvalueOrd = pvalue[order(pvalue)[1:numPlot]]

    if(plot==TRUE) {
        par(las=2)
        par(mar=c(9,4,1,0.5))
        b <- barplot(-log10(pvalueOrd),main=main, ylim=c(0, -log10(min(pvalueOrd))*1.1), ylab=expression(-log[10](pvalue)), xaxt = "n")
        text(x=b,y=-log10(pvalueOrd),labels=signif(pvalueOrd,2),pos=3)
        text(x=b, y=log10(pvalueOrd[1])*0.02, names(pvalueOrd), xpd=TRUE, srt=45, pos=2)
    }

    pvalue
}

# Read pair of .cov files from bismark_methylation_extractor to build table methylSig expects
# fileList is a list of files where each entry in the list consists of a vector pair of files
# *bismark.cov is first and *cytosine.cov is second.
# Called by readBismarkData
readBismarkOutputSingleFile = function(fileIndex, fileList, minCount, maxCount, destranded, filterSNPs, quiet=FALSE) {
    # Read files and minimize memory footprint with colClasses
    message(sprintf('Reading %s',fileList[[fileIndex]][1]))
    cov = read.table(fileList[[fileIndex]][1], sep='\t', header=F,
        col.names=c('chr','start','end','perc_meth','numCs','numTs'),
        colClasses=c('character','numeric','NULL','NULL','numeric','numeric'),stringsAsFactors=F)

    message(sprintf('Reading %s',fileList[[fileIndex]][2]))
    cyt = read.table(fileList[[fileIndex]][2], sep='\t', header=F,
        col.names=c('chr','pos','strand','numCs','numTs','C_context','tri_context'),
        colClasses=c('character','numeric','character','NULL','NULL','NULL','NULL'), stringsAsFactors=F)

    # Create the chromBase column which matching will be done on
    cov$chromBase = paste(cov$chr, cov$start, sep='.')
    cyt$chromBase = paste(cyt$chr, cyt$pos, sep='.')

    # Extract strand from the cytosine report
    # This is the only purpose of the cytosine report
    cov$strand = cyt[match(cov$chromBase,cyt$chromBase), 'strand']

    # Remove cyt to minimize memory footprint
    rm(cyt)

    # Indicate if some sites in *bismark.cov are not present in *cytosine.cov
    # This probably shouldn't happen, but it would mean the strand is NA.
    if(length(which(is.na(cov$strand))) > 0) {message('WARNING! Positions in coverage file not found in cytosine report.')}

    cov$coverage = cov$numCs + cov$numTs

    # If the sum of frequencies of C's and T's is less than 95, discard sites
    # by setting coverage to 0, and report for each pair of files the
    # proportion of sites that are removed.
    invalidList = which( ((cov$numCs + cov$numTs) / cov$coverage) < 0.95)
    if(length(invalidList) > 0) {
        cat("(", fileIndex,"/", NROW(fileList), ") ", "Frequency Invalid List: ", NROW(invalidList), "/", NROW(cov), "=", signif(NROW(invalidList)/NROW(cov),3), "\n", sep="")
        cov$coverage[invalidList] = 0
    }

    # Set coverage of sites exceeding maxCount or not exceeding minCount to 0
    invalidList = which(cov$coverage > maxCount | cov$coverage < minCount)
    cat("(", fileIndex,"/", NROW(fileList), ") ", "Count Invalid List: ", NROW(invalidList), "/", NROW(cov), "=", signif(NROW(invalidList)/NROW(cov),3), "\n", sep="")
    cov$coverage[invalidList] = 0

    # In case of destranded, shift reverse strand CpG sites to match forward strand
    if(destranded == TRUE) {
        invalidList = which(cov$strand == "-" | cov$strand == "R")
        cov$start[invalidList] = cov$start[invalidList] - 1
        # The end column is ignored
        # cov$end[invalidList] = cov$end[invalidList] - 1
    }

    if(filterSNPs) {
        data('CT_SNPs_hg19',envir=environment())
        cov_gr = GRanges(seqnames=cov$chr, ranges=IRanges(start=cov$start, end=cov$start))

        overlaps = findOverlaps(cov_gr, CT_SNPs_hg19)
        invalidList = overlaps@queryHits

        cat("(", fileIndex,"/", NROW(fileList), ") ", "SNP Invalid List: ", NROW(invalidList), "/", NROW(cov), "=", signif(NROW(invalidList)/NROW(cov),3), "\n", sep="")

        cov$coverage[invalidList] = 0
    }

    # Pull out and rename relevant columns
    final = cov[,c('chromBase','chr','start','strand','coverage','numCs','numTs')]
    # Matching column names from methylSigReadData
    colnames(final) = c('id','chr','start','strand','coverage','numCs','numTs')

    # Make strand column a factor.
    final$strand = as.factor(final$strand)
    ## When only F observed in strand column, as.factor convert F to FALSE
    levels(final$strand) = list("+"="F","-"="R","*"="*", "+"="FALSE")

    final = final[final$coverage > 0,]

    return(final)
}

#' Read output from bismark_methylation_extractor to make a 'methylSigData' object.
#'
#' This function takes the coverage and cytosine report files from the \code{bismark_methylation_extractor} (options \code{--bedGraph} and \code{--cytosine_report}) and constructs a table conforming to that expected by the 'methylSigData' class.
#'
#' @param bismarkCovFiles Vector of coverage files (.cov as of Bismark v0.13.0) from \code{bismark_methylation_extractor} with \code{--bedGraph} flag. Sample order should match that of \code{cytosineCovFiles}.
#' @param cytosineCovFiles Vector of cytosine report files from \code{bismark_methylation_extractor} with \code{--cytosine_report} flag. Sample order should match that of \code{bismarkCovFiles}.
#' @param sample.ids Vector of sample ids.
#' @param assembly Character string indicating the genome assembly, such as "hg18", "hg19", "mm9", or "mm10".
#' @param pipeline Character string indicating the pipepline name that generated the data, for example, "bismark".
#' @param context Methylation context string such "CpG","CpH",or "CHH".
#' @param resolution A string indicating whether the input data are  base-pair or regional resolution. Either "base" or "region" is allowed. Default is "base".
#' @param treatment A numeric vector contatining numbers to distinguish the group classification.
#' @param destranded A logical value indicating whether to destrand the reverse to forward strand. If TRUE, the reads from both will be combined. Default is TRUE.
#' @param maxCount A number indicating the maximum coverage count to be included.
#' @param minCount A number indicating the minimum coverage count to be included.
#' @param filterSNPs A logical value indicating whether or not to filter out C > T SNPs based on the 1000 Genomes Project.
#' @param num.cores Number of cores to be used in reading files. Default is 1.
#' @param quiet A logical value. If quiet=TRUE, then this function does not show progress information. Default is FALSE.
#'
#' @return A 'methylSigData' object.
#'
#' @seealso \code{\link{methylSigCalc}}
#'
#' @export
readBismarkData = function(bismarkCovFiles, cytosineCovFiles,
            sample.ids, assembly=NA, pipeline=NA, context=NA,resolution="base",treatment,
            destranded=TRUE, maxCount=500, minCount=10, filterSNPs=FALSE, num.cores=1, quiet=FALSE) {

    # Checks to verify *bismark.cov files match *cytosine.cov files
    if(length(bismarkCovFiles) != length(cytosineCovFiles)){
        stop('The number of *bismark.cov files does not match the number of *cytosine.cov files! Each sample should have these files from bismark_methylation_extractor.')
    }

    fileList = lapply(1:length(bismarkCovFiles),function(i){c(bismarkCovFiles[i],cytosineCovFiles[i])})

    n.files = length(fileList)

    if(num.cores > 1) {
        chrList <- mclapply(1:n.files, readBismarkOutputSingleFile, fileList, minCount=minCount, maxCount=maxCount, destranded, filterSNPs, quiet, mc.cores=num.cores)
    } else {
        chrList = list()
        for(i in 1:n.files) chrList[[i]] =  readBismarkOutputSingleFile(i, fileList, minCount=minCount, maxCount=maxCount, destranded, filterSNPs, quiet)
    }

    MAXBASE = 0
    uniqueChr = NULL
    for(fileIndex in 1:n.files) {
         uniqueChr = c(uniqueChr, chrList[[fileIndex]]$chr)
         MAXBASE = max(MAXBASE, max(chrList[[fileIndex]]$start))
    }

    uniqueChr = unique(uniqueChr)
    uniqueChr = uniqueChr[order(uniqueChr)]

    # This is kind of a hash?
    MAXBASE = 10^{ceiling(log10(MAXBASE + 1))}
    uniqueLoc = NULL
    for(fileIndex in 1:n.files) {
        chrList[[fileIndex]]$chr = factor(chrList[[fileIndex]]$chr, levels=uniqueChr)
        uniqueLoc = unique(c(uniqueLoc, as.numeric(chrList[[fileIndex]]$chr) * MAXBASE+chrList[[fileIndex]]$start))
    }
    uniqueLoc = uniqueLoc[order(uniqueLoc)]

    # This dictates the sites in common across the samples?
    sizeRet = NROW(uniqueLoc)

    coverage = numCs = numTs = matrix(0, nrow=sizeRet, ncol=n.files)
    strand = factor(rep(NA, sizeRet), levels=levels(chrList[[1]]$strand))

    for(fileIndex in 1:n.files) {
        if(quiet == FALSE) {
            cat("(",fileIndex,")",sep="")
            if(fileIndex %% 10 == 0) cat("\n")
        }
        location =  findInterval( (as.numeric(chrList[[fileIndex]]$chr) * MAXBASE+chrList[[fileIndex]]$start) , uniqueLoc)
        if(destranded == FALSE) {
            strand[location] = chrList[[fileIndex]]$strand
            coverage[location,fileIndex] = chrList[[fileIndex]]$coverage
            numCs[location,fileIndex] = chrList[[fileIndex]]$numCs
            numTs[location,fileIndex] = chrList[[fileIndex]]$numTs
        } else {
            settingList = is.na(strand[location])
            strand[location][settingList] = chrList[[fileIndex]]$strand[settingList]
            settingList = !settingList & (strand[location] != chrList[[fileIndex]]$strand)
            strand[location][settingList] = "*"

            forward = (chrList[[fileIndex]]$strand == "+")
            coverage[location[forward],fileIndex] = chrList[[fileIndex]]$coverage[forward]
            numCs[location[forward],fileIndex] = chrList[[fileIndex]]$numCs[forward]
            numTs[location[forward],fileIndex] = chrList[[fileIndex]]$numTs[forward]
            reverse = (chrList[[fileIndex]]$strand == "-")
            coverage[location[reverse],fileIndex] = coverage[location[reverse],fileIndex]+chrList[[fileIndex]]$coverage[reverse]
            numCs[location[reverse],fileIndex] = numCs[location[reverse],fileIndex] + chrList[[fileIndex]]$numCs[reverse]
            numTs[location[reverse],fileIndex] = numTs[location[reverse],fileIndex] + chrList[[fileIndex]]$numTs[reverse]

        }
    }

    if(quiet == FALSE) cat("\n")

#    strand = as.factor(strand)
#    levels(strand) = list("+"="1","-"="2","*"="3")

    # Convert fileList to a character vector
    fileList = unlist(fileList)

    options = paste("maxCount=", maxCount, " & minCount=", minCount, " & filterSNPs=", filterSNPs, sep="")
    if(!is.na(assembly)) options = paste(options, " & assembly=", assembly, sep="")
    if(!is.na(context))  options = paste(options, " & context=",  context,  sep="")
    if(!is.na(pipeline)) options = paste(options, " & pipeline=", pipeline, sep="")

    numTs[coverage==0] = NA
    numCs[coverage==0] = NA
    coverage[coverage==0] = NA
    methylSig.newData(data.ids=uniqueLoc, data.chr=as.factor(uniqueChr[as.integer(uniqueLoc/MAXBASE)]),
                      data.start=uniqueLoc%%MAXBASE, data.end=uniqueLoc%%MAXBASE,
                      data.strand=strand, data.coverage = coverage, data.numTs = numTs, data.numCs = numCs,
                      sample.ids=sample.ids, treatment=treatment, destranded=destranded,
                      resolution=resolution, sample.filenames=fileList,options=options)
}

# Called by methylSigReadData
methylSigReadDataSingleFile <- function(fileIndex, fileList, header, minCount, maxCount, destranded, filterSNPs, quiet=FALSE) {
    if(quiet==FALSE) cat("Reading file (", fileIndex, "/", NROW(fileList),") -- ", fileList[[fileIndex]], "\n", sep="")
    chr = read.table(fileList[[fileIndex]], header=header, stringsAsFactors=FALSE)
    #### order for base ####
    ##chr = chr[order(chr$base),]
    ####
    names(chr) <- c("id",  "chr", "start", "strand", "coverage", "numCs", "numTs")

    # At this point, numCs and numTs are actually freqC and freqT from methylKit.
    invalidList = which(chr$numCs + chr$numTs < 95)
    if(length(invalidList) > 0) {
        cat("(", fileIndex,"/", NROW(fileList), ") ", "Frequency Invalid List: ", NROW(invalidList), "/", NROW(chr), "=", signif(NROW(invalidList)/NROW(chr),3), "\n", sep="")
        chr$coverage[invalidList] = 0
    }

    # Now numCs and numTs have frequencies replaced by counts
    chr$numCs<- round(chr$numCs * chr$coverage / 100)
    chr$numTs<- round(chr$numTs * chr$coverage / 100)

    invalidList = which(chr$coverage > maxCount | chr$coverage < minCount)
    cat("(", fileIndex,"/", NROW(fileList), ") ", "Count Invalid List: ", NROW(invalidList), "/", NROW(chr), "=", signif(NROW(invalidList)/NROW(chr),3), "\n", sep="")
    chr$coverage[invalidList] = 0

    if(destranded == TRUE) {
        invalidList = which(chr$strand == "-" | chr$strand == "R")
        chr$start[invalidList] = chr$start[invalidList] - 1
        # The files coming from methylKit don't have end columns
        # chr$end[invalidList] = chr$end[invalidList] - 1
    }

    if(filterSNPs) {
        data('CT_SNPs_hg19',envir=environment())
        chr_gr = GRanges(seqnames=chr$chr, ranges=IRanges(start=chr$start, end=chr$start))

        overlaps = findOverlaps(chr_gr, CT_SNPs_hg19)
        invalidList = overlaps@queryHits

        cat("(", fileIndex,"/", NROW(fileList), ") ", "SNP Invalid List: ", NROW(invalidList), "/", NROW(chr), "=", signif(NROW(invalidList)/NROW(chr),3), "\n", sep="")

        chr$coverage[invalidList] = 0
    }

    #chr$chr = as.factor(chr$chr)
    chr$strand = as.factor(chr$strand)
    ## When only F observed in strand column, as.factor convert F to FALSE
    levels(chr$strand) = list("+"="F","-"="R","*"="*", "+"="FALSE")

    chr[chr$coverage>0,2:7]
}

#' Read methylation score files to make a 'methylSigData' object.
#'
#' This function reads methylation score files (having columns chrBase, chr, base, strand, coverage, freqC, freqT) to make a \code{\link{methylSigData-class}} object that can be used in differential methylation analysis.
#'
#' @param fileList Vector of files to be read. The methylKit package can be used to generate the appropriate tables from sorted SAM output from Bismark.
#' @param sample.ids Vector of sample ids.
#' @param assembly Character string indicating the genome assembly, such as "hg18", "hg19", "mm9", or "mm10".
#' @param pipeline Character string indicating the pipepline name that generated the data, for example, "bismark".
#' @param header A logical value indicating whether the score files have header or not.  Default is TRUE.
#' @param context Methylation context string such "CpG","CpH",or "CHH".
#' @param resolution A string indicating whether the input data are  base-pair or regional resolution. Either "base" or "region" is allowed. Default is "base".
#' @param treatment A numeric vector contatining numbers to distinguish the group classification.
#' @param destranded A logical value indicating whether to destrand the reverse to forward strand. If TRUE, the reads from both will be combined. Default is TRUE.
#' @param maxCount A number indicating the maximum coverage count to be included.
#' @param minCount A number indicating the minimum coverage count to be included.
#' @param filterSNPs A logical value indicating whether or not to filter out C > T SNPs based on the 1000 Genomes Project.
#' @param num.cores Number of cores to be used in reading files. Default is 1.
#' @param quiet A logical value. If quiet=TRUE, then this function does not show progress information. Default is FALSE.
#'
#' @return A \code{\link{methylSigData-class}} object.
#'
#' @seealso \code{\link{methylSigCalc}}
#'
#' @examples
#' fileList = c( system.file("extdata", "AML_1.txt", package = "methylSig"),
#'               system.file("extdata", "AML_2.txt", package = "methylSig"),
#'               system.file("extdata", "AML_3.txt", package = "methylSig"),
#'               system.file("extdata", "AML_4.txt", package = "methylSig"),
#'               system.file("extdata", "NBM_1.txt", package = "methylSig"),
#'               system.file("extdata", "NBM_2.txt", package = "methylSig"),
#'               system.file("extdata", "NBM_3.txt", package = "methylSig"),
#'               system.file("extdata", "NBM_4.txt", package = "methylSig"))
#' sample.id = c("AML1","AML2","AML3","AML4","NBM1","NBM2","NBM3","NBM4")
#'
#' treatment = c(1,1,1,1,0,0,0,0)
#'
#' meth <- methylSigReadData(fileList, sample.ids = sample.id,
#'             assembly = "hg18", treatment = treatment,
#'             context = "CpG", destranded=TRUE)
#'
#' @export
methylSigReadData = function(fileList,
            sample.ids, assembly=NA, pipeline=NA, header=TRUE, context=NA,resolution="base",treatment,
            destranded=TRUE, maxCount=500, minCount=10, filterSNPs=FALSE, num.cores=1, quiet=FALSE) {

    n.files = NROW(fileList)

    if(num.cores > 1) {
        chrList <- mclapply(1:n.files, methylSigReadDataSingleFile, fileList, header = header, minCount=minCount, maxCount=maxCount, destranded, filterSNPs, quiet, mc.cores=num.cores)
    } else {
        chrList = list()
        for(i in 1:n.files) chrList[[i]] =  methylSigReadDataSingleFile(i, fileList, header=header, minCount=minCount, maxCount=maxCount, destranded, filterSNPs, quiet)
    }


    MAXBASE = 0
    uniqueChr = NULL
    for(fileIndex in 1:n.files) {
         uniqueChr = c(uniqueChr, chrList[[fileIndex]]$chr)
         MAXBASE = max(MAXBASE, max(chrList[[fileIndex]]$start))
    }

    uniqueChr = unique(uniqueChr)
    uniqueChr = uniqueChr[order(uniqueChr)]

    MAXBASE = 10^{ceiling(log10(MAXBASE + 1))}
    uniqueLoc = NULL
    for(fileIndex in 1:n.files) {
        chrList[[fileIndex]]$chr = factor(chrList[[fileIndex]]$chr, levels=uniqueChr)
        uniqueLoc = unique(c(uniqueLoc, as.numeric(chrList[[fileIndex]]$chr) * MAXBASE+chrList[[fileIndex]]$start))
    }

    uniqueLoc = uniqueLoc[order(uniqueLoc)]

    sizeRet = NROW(uniqueLoc)

    coverage = numCs = numTs = matrix(0, nrow=sizeRet, ncol=n.files)
    strand = factor(rep(NA, sizeRet), levels=levels(chrList[[1]]$strand))

    for(fileIndex in 1:n.files) {
        if(quiet == FALSE) {
            cat("(",fileIndex,")",sep="")
            if(fileIndex %% 10 == 0) cat("\n")
        }
        location =  findInterval( (as.numeric(chrList[[fileIndex]]$chr) * MAXBASE+chrList[[fileIndex]]$start) , uniqueLoc)
        if(destranded == FALSE) {
            strand[location] = chrList[[fileIndex]]$strand
            coverage[location,fileIndex] = chrList[[fileIndex]]$coverage
            numCs[location,fileIndex] = chrList[[fileIndex]]$numCs
            numTs[location,fileIndex] = chrList[[fileIndex]]$numTs
        } else {
            settingList = is.na(strand[location])
            strand[location][settingList] = chrList[[fileIndex]]$strand[settingList]
            settingList = !settingList & (strand[location] != chrList[[fileIndex]]$strand)
            strand[location][settingList] = "*"

            forward = (chrList[[fileIndex]]$strand == "+")
            coverage[location[forward],fileIndex] = chrList[[fileIndex]]$coverage[forward]
            numCs[location[forward],fileIndex] = chrList[[fileIndex]]$numCs[forward]
            numTs[location[forward],fileIndex] = chrList[[fileIndex]]$numTs[forward]
            reverse = (chrList[[fileIndex]]$strand == "-")
            coverage[location[reverse],fileIndex] = coverage[location[reverse],fileIndex]+chrList[[fileIndex]]$coverage[reverse]
            numCs[location[reverse],fileIndex] = numCs[location[reverse],fileIndex] + chrList[[fileIndex]]$numCs[reverse]
            numTs[location[reverse],fileIndex] = numTs[location[reverse],fileIndex] + chrList[[fileIndex]]$numTs[reverse]

        }
    }

    if(quiet == FALSE) cat("\n")

#    strand = as.factor(strand)
#    levels(strand) = list("+"="1","-"="2","*"="3")

    options = paste("maxCount=", maxCount, " & minCount=", minCount, " & filterSNPs=", filterSNPs, sep="")
    if(!is.na(assembly)) options = paste(options, " & assembly=", assembly, sep="")
    if(!is.na(context))  options = paste(options, " & context=",  context,  sep="")
    if(!is.na(pipeline)) options = paste(options, " & pipeline=", pipeline, sep="")

    numTs[coverage==0] = NA
    numCs[coverage==0] = NA
    coverage[coverage==0] = NA
    methylSig.newData(data.ids=uniqueLoc, data.chr=as.factor(uniqueChr[as.integer(uniqueLoc/MAXBASE)]),
                      data.start=uniqueLoc%%MAXBASE, data.end=uniqueLoc%%MAXBASE,
                      data.strand=strand, data.coverage = coverage, data.numTs = numTs, data.numCs = numCs,
                      sample.ids=sample.ids, treatment=treatment, destranded=destranded,
                      resolution=resolution, sample.filenames=fileList,options=options)
}

# Not called
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
#' @param chr A single string indicating chromosome from `methylSigData' or `methylSigDiff' objects.
#' @param start A vector of CpG site loci from `methylSigData' or `methylSigDiff' objects.
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
