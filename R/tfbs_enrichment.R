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

#' Perform transcription factor enrichment test among differentially methylated cytosines or regions
#'
#' This function tests for enriched transcription binding sites among differentially methylated sites or regions using a binomial test. It also has an option to generate a bar plot to show the significance levels for the enriched transcription factor binding sites.
#'
#' Likelihood ratio test is used based on the binomial distribution.
#'
#' @param myDiff \code{methylSigDiff-class} object that contains all CpG sites that are tested for differential methylation.
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
