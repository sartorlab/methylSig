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

#' Annotation function using RefGene models
#'
#' This function annotates the CpG sites to promoter, coding DNA sequence (CDS), exons, 3' untranslated region, 5' untranslated region, nocoding RNA region, or integenic regions using the RefGene gene models.
#'
#'   This function generates an object to annotate refGene information. The same genome assembly should be used for `myDiff' and `refGeneInfo'. The promoter region is defined as 1,000bp upstream from a transcription starting site (TSS) to the TSS.
#'
#' @param refGeneInfo RefGene information object from function \code{\link{getRefgeneInfo}}.
#' @param myDiff \code{methylSigDiff-class} object to be used for annotation.
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

#' CpG island annotation function
#'
#' This function annotates the CpG sites or regions to CpG island, CpG shore, CpG shelf or interCGI regions.
#'
#'   This function generates an object to annotate CpG island information. The same genome assembly should be used for `myDiff' and `cpgInfo'. CpG shores are defined as the region outside CpG islands but within 2,000 bp of any CpG island. CpG shelves are defined as the region within 2,000 bp away from a CpG shore. When regions overlap, the priority is CpG island, then CpG shore. Regions >4,000 bp from a CpG island are deemed interCGI.
#'
#' @param cpgInfo A CpG island information object.
#' @param myDiff A \code{methylSigDiff-class} object to be used for annotation.
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
