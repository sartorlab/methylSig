
#####################################################
##
##  Dec. 23th 2013
##
#####################################################

##### May 20 2014 ####################
write.methylSigDiff <- function(object, ...) {
    printRange = 1:NROW(object@data.ids)
    printData = data.frame(chr=object@data.chr[printRange], start=object@data.start[printRange],
                           end=object@data.end[printRange], strand=object@data.strand[printRange],
                           object@results[printRange,,drop=FALSE])
    write.table(printData, ...)
#    cat("sample.ids:",object@sample.ids,"\n", sep=" ")
#    cat("treatment:", object@treatment,"\n", sep=" ")
#    cat("destranded:",object@destranded,"\n", sep=" ")
#    cat("resolution:", object@resolution,"\n", sep=" ")
#    cat("options:", object@options,"\n", sep=" ")
}


##############                ################
############## Binomial Model ################
##############(Dec. 23, 2013) ################

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


############################################
####           Dec. 24, 2013           #####
############################################

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
                      resolution="region",sample.filenames=meth@sample.filenames,options=meth@options)
}

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
                      resolution="TF",sample.filenames=meth@sample.filenames,options="")
}


######## Functions ########
#                         #
#    Row: Samples         #
#    Column: Locations    #
###########################

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
methylSig_weightFunc <- function(u) (1-u^2)^3

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
                              resolution=meth@resolution, options=options, data.options = meth@options)
}


#############################
###  loc.range= from, to  ###
#############################

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

getTFBSCountByChrom <- function(tfbsInfo, startEnd, listToCount) {
###
    maxOverLaps = 200
    cpgId = findInterval(listToCount,tfbsInfo[[2]]$chromStart[startEnd])
    cpgIdAll = pmax(0,rep(cpgId,each=maxOverLaps)+rep(-(maxOverLaps-1):0,NROW(cpgId)))
    whichValid = which(rep(listToCount,each=maxOverLaps)<=c(0,tfbsInfo[[2]]$chromEnd[startEnd])[cpgIdAll+1])

#### In case of that not all slots have data   
    table(c(1:NROW(levels(tfbsInfo[[2]]$name)),tfbsInfo[[2]]$name[startEnd][cpgIdAll[whichValid]]))-1
}

is.TFBSByChrom <- function(tfbsInfo, startEnd, listToCount) {
###
    maxOverLaps = 200
    cpgId = findInterval(listToCount,tfbsInfo[[2]]$chromStart[startEnd])
    cpgIdAll = pmax(0,rep(cpgId,each=maxOverLaps)+rep(-(maxOverLaps-1):0,NROW(cpgId)))
    whichValid = rep(listToCount,each=maxOverLaps)<=c(0,tfbsInfo[[2]]$chromEnd[startEnd])[cpgIdAll+1]

#### In case of that not all slots have data   
    colSums(matrix(whichValid, nrow=maxOverLaps))>0
}

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

refGeneAnnotationPlot <- function(listFrom, main="Plot", priority=c("cds", "promoter","noncoding", "5'utr", "3'utr")) {
    countInfo= rep(NA, NCOL(listFrom))
#####
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


########################################
###     Revised Dec. 23, 2012      #####
########################################

methylSigReadDataSingleFile <- function(fileIndex, fileList, header, minCount, maxCount, destranded, quiet=FALSE) {
    if(quiet==FALSE) cat("Reading file (", fileIndex, "/", NROW(fileList),") -- ", fileList[[fileIndex]], "\n", sep="")
    chr = read.table(fileList[[fileIndex]], header=header, stringsAsFactors=FALSE)
    #### order for base ####
    ##chr = chr[order(chr$base),]
#### 
    names(chr) <- c("id",  "chr", "start", "strand", "coverage", "numCs", "numTs")

    invalidList = which(chr$numCs + chr$numTs < 95) 
    if(length(invalidList) > 0) {
        cat("(", fileIndex,"/", NROW(fileList), ") ", "Invalid List: ", NROW(invalidList), "/", NROW(chr), "=", signif(NROW(invalidList)/NROW(chr),3), "\n", sep="")
        chr$coverage[invalidList] = 0
    }

    chr$numCs<- round(chr$numCs * chr$coverage / 100)
    chr$numTs<- round(chr$numTs * chr$coverage / 100)
    size = NROW(chr$base)

    invalidList = (chr$coverage > maxCount | chr$coverage < minCount)
    chr$coverage[invalidList] = 0

    if(destranded == TRUE) {
        invalidList = which(chr$strand == "-" | chr$strand == "R")
        chr$start[invalidList] = chr$start[invalidList] - 1
        chr$end[invalidList] = chr$end[invalidList] - 1
    } 

    #chr$chr = as.factor(chr$chr)
    chr$strand = as.factor(chr$strand)
    ## When only F observed in strand column, as.factor convert F to FALSE
    levels(chr$strand) = list("+"="F","-"="R","*"="*", "+"="FALSE")

    chr[chr$coverage>0,2:7]
}

### Revised Dec. 23, 2012
###
methylSigReadData = function(fileList,
            sample.ids, assembly=NA, pipeline=NA, header=TRUE, context=NA,resolution="base",treatment,
            destranded=TRUE, maxCount=500, minCount=10, num.cores=1, quiet=FALSE) {

    n.files = NROW(fileList)

    if(num.cores > 1) {
        chrList <- mclapply(1:n.files, methylSigReadDataSingleFile, fileList, header = header, minCount=minCount, maxCount=maxCount, destranded, quiet, mc.cores=num.cores)
    } else {
        chrList = list()
        for(i in 1:n.files) chrList[[i]] =  methylSigReadDataSingleFile(i, fileList, header=header, minCount=minCount, maxCount=maxCount, destranded, quiet)
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
        uniqueLoc = unique(c(uniqueLoc,as.numeric(chrList[[fileIndex]]$chr)*MAXBASE+chrList[[fileIndex]]$start))
    }

    uniqueLoc = uniqueLoc[order(uniqueLoc)]   

    sizeRet = NROW(uniqueLoc)

    coverage=numCs=numTs = matrix(0, nrow=sizeRet, ncol=n.files)
    strand = factor(rep(NA, sizeRet), levels=levels(chrList[[1]]$strand))

    for(fileIndex in 1:n.files) {
        if(quiet == FALSE) {
            cat("(",fileIndex,")",sep="")
            if(fileIndex %% 10 == 0) cat("\n")
        }
        location =  findInterval((as.numeric(chrList[[fileIndex]]$chr)*MAXBASE+chrList[[fileIndex]]$start),uniqueLoc)
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

    options = paste("maxCount=", maxCount, " & minCount=", minCount, sep="")
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
                      resolution=resolution,sample.filenames=fileList,options=options)
} 

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

