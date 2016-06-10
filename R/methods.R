#' This class is for generating 'methylSigDiff' object.
#'
#' Objects can be created by calls of the form \code{new("methylSigDiff", ...)}.
#'
#' @slot data.ids A numeric vector
#' @slot data.chr A factor vector
#' @slot data.strand A factor vector
#' @slot data.start A numeric vector
#' @slot data.end A numeric vector
#' @slot results A matrix
#' @slot treatment A numeric vector
#' @slot sample.ids A character vector
#' @slot sample.filenames A character vector
#' @slot destranded A logical value
#' @slot resolution A character value
#' @slot options A character vector
#' @slot data.options A character vector
#'
#' @keywords classes internal
#' @aliases methylSigDiff-class methylSigDiff
#'
#' @export
setClass("methylSigDiff", representation(data.ids="numeric", data.chr="factor", data.start="numeric",
                                         data.end="numeric", data.strand="factor", results="matrix",
                                         sample.ids = "character", treatment="numeric", sample.filenames = "character",
                                         destranded="logical", resolution="character",
                                         options="character", data.options="character"))

methylSig.newDiff <- function(data.ids, data.chr, data.start,data.end,
    data.strand, results, sample.ids="", sample.filenames="", treatment,
    destranded=TRUE, resolution="", options="", data.options="") {
    ret = new("methylSigDiff")
    ret@data.ids=data.ids
    ret@data.chr = data.chr
    ret@data.strand=data.strand
    ret@data.start=data.start
    ret@data.end=data.end
    ret@results = results
    ret@treatment = treatment
    ret@sample.ids=sample.ids
    ret@sample.filenames = sample.filenames
    ret@destranded = destranded
    ret@resolution=resolution
    ret@options=options
    ret@data.options=data.options

    ret
}

methylSig.subDiff <- function(meth,i) {
    ret = new("methylSigDiff")
    ret@data.ids=meth@data.ids[i,drop=FALSE]
    ret@data.chr = meth@data.chr[i,drop=FALSE]
    ret@data.strand=meth@data.strand[i,drop=FALSE]
    ret@data.start=meth@data.start[i,drop=FALSE]
    ret@data.end=meth@data.end[i,drop=FALSE]
    ret@results = meth@results[i,,drop=FALSE]
    ret@treatment = meth@treatment
    ret@sample.ids=meth@sample.ids
    ret@sample.filenames = meth@sample.filenames
    ret@destranded = meth@destranded
    ret@resolution=meth@resolution
    ret@options=meth@options
    ret@data.options=meth@data.options

    ret
}

methlSig.subDiffElement<-function(meth,i,j) {
    if(j=="chr") {
        meth@data.chr[i]
    } else if(j=="start") {
        meth@data.start[i]
    } else if(j=="end") {
        meth@data.end[i]
    } else if(j=="strand") {
        meth@data.strand[i]
    } else {
        meth@results[i,j]
    }
}

#### sub-setting ####
### meth[1:10,1:9] : obtain first 10 lines from first 9 samples
###

#' Subset methylSigDiff object
#'
#' @param i Start position
#' @param j End position
#'
#' @rdname methylSigDiff-class
#'
#' @keywords methods internal
setMethod("[", signature(x = "methylSigDiff", i = "ANY", j = "character"),
    function(x, i, j) {
        if(missing(i)) i=1:NROW(x@data.chr)
        methlSig.subDiffElement(x,i,j)
    }
)

#' Subset methylSigDiff object
#'
#' @param i Start position
#' @param j End position
#'
#' @rdname methylSigDiff-class
#'
#' @keywords methods internal
setMethod("[", signature(x = "methylSigDiff", i = "ANY", j = "missing"),
    function(x, i) {
       if(missing(i)) {
           x
       } else methylSig.subDiff(x,i)
    }
)

#' This method is for printing a methylSigDiff object.
#'
#' @param object A methylSigDiff object.
#'
#' @rdname methylSigDiff-class
#'
#' @examples
#' show("myDiff")
#'
#' @keywords methods internal
setMethod("show", "methylSigDiff",
  function(object) {
#    pvalue=0,qvalue=0,meth.diff=0,logLikRatio=NA, phi=0, df=0, muT=0, muC=0)

    cat("methylSigDiff object with", format(NROW(object@data.ids), big.mark=",",small.interval=3), "rows\n")
    cat("--------------------------\n")
    printRange = 1:min(10,NROW(object@data.ids))
    printData = data.frame(chr=object@data.chr[printRange], start=object@data.start[printRange],
                           end=object@data.end[printRange], strand=object@data.strand[printRange],
                           object@results[printRange,,drop=FALSE])
    print(printData)
    cat("--------------------------\n")
    cat("sample.ids:",object@sample.ids,"\n", sep=" ")
    cat("treatment:", object@treatment,"\n", sep=" ")
    cat("destranded:",object@destranded,"\n", sep=" ")
    cat("resolution:", object@resolution,"\n", sep=" ")
    cat("options:", object@options,"\n", sep=" ")
  }
)

#' This class is for generating 'methylSigData' object.
#'
#' Objects can be created by calls of the form \code{new("methylSigData", ...)}.
#'
#' @slot data.ids A numeric vector
#' @slot data.chr A factor vector
#' @slot data.strand A factor vector
#' @slot data.start A numeric vector
#' @slot data.end A numeric vector
#' @slot data.coverage A matrix
#' @slot data.numCs A matrix
#' @slot data.numTs A matrix
#' @slot treatment A numeric vector
#' @slot sample.ids A character vector
#' @slot sample.filenames A character vector
#' @slot destranded A logical value
#' @slot resolution A character value
#' @slot options A character vector
#'
#' @examples
#' showClass("methylSigData")
#'
#' @aliases methylSigData-class methylSigData
#' @rdname methylSigData-class
#'
#' @keywords classes internal
#'
#' @export
setClass("methylSigData", representation(data.ids="numeric", data.chr="factor", data.start="numeric",
                                         data.end="numeric", data.strand="factor", data.coverage="matrix",
                                         data.numCs = "matrix", data.numTs = "matrix",
                                         sample.ids = "character", sample.filenames = "character", treatment="numeric",
                                         destranded="logical", resolution="character",
                                         options="character"))

methylSig.newData <- function(data.ids, data.chr, data.start,data.end,
    data.strand, data.coverage, data.numCs, data.numTs, sample.ids="",
    sample.filenames="", treatment, destranded=TRUE, resolution="",
    options="") {
    ret = new("methylSigData")
    ret@data.ids=data.ids
    ret@data.chr = data.chr
    ret@data.strand=data.strand
    ret@data.start=data.start
    ret@data.end=data.end
    ret@data.coverage=data.coverage
    ret@data.numCs = data.numCs
    ret@data.numTs = data.numTs
    ret@treatment = treatment
    ret@sample.ids=sample.ids
    ret@sample.filenames = sample.filenames
    ret@destranded = destranded
    ret@resolution=resolution
    ret@options=options

    ret
}

methylSig.subData <- function(meth,i,j, d) {
    ret = new("methylSigData")

    if(d) {
        nonZeroRow = i[rowSums(meth@data.coverage[i,j,drop=FALSE] > 0, na.rm=TRUE) > 0]

#        nonZeroRow = which(rowSums(ret@data.coverage > 0, na.rm=TRUE) > 0)
        ret@data.ids=meth@data.ids[nonZeroRow,drop=FALSE]
        ret@data.chr = meth@data.chr[nonZeroRow,drop=FALSE]
        ret@data.strand=meth@data.strand[nonZeroRow,drop=FALSE]
        ret@data.start=meth@data.start[nonZeroRow,drop=FALSE]
        ret@data.end=meth@data.end[nonZeroRow,drop=FALSE]
        ret@data.coverage=meth@data.coverage[nonZeroRow,j,drop=FALSE]
        ret@data.numCs=meth@data.numCs[nonZeroRow,j,drop=FALSE]
        ret@data.numTs=meth@data.numTs[nonZeroRow,j,drop=FALSE]
    } else {
        ret@data.ids=meth@data.ids[i,drop=FALSE]
        ret@data.chr = meth@data.chr[i,drop=FALSE]
        ret@data.strand=meth@data.strand[i,drop=FALSE]
        ret@data.start=meth@data.start[i,drop=FALSE]
        ret@data.end=meth@data.end[i,drop=FALSE]
        ret@data.coverage=meth@data.coverage[i,j,drop=FALSE]
        ret@data.numCs = meth@data.numCs[i,j,drop=FALSE]
        ret@data.numTs = meth@data.numTs[i,j,drop=FALSE]
    }

    ret@sample.filenames = meth@sample.filenames[j,drop=FALSE]
    ret@treatment = meth@treatment[j,drop=FALSE]
    ret@sample.ids=meth@sample.ids[j,drop=FALSE]

    ret@destranded = meth@destranded
    ret@resolution=meth@resolution
    ret@options=meth@options

    ret
}

methylSig.subDataElement<-function(meth,i,j) {
    if(j=="chr") {
        meth@data.chr[i]
    } else if(j=="start") {
        meth@data.start[i]
    } else if(j=="end") {
        meth@data.end[i]
    } else if(j=="strand") {
        meth@data.strand[i]
    } else if(j=="coverage") {
        meth@data.coverage
    } else if(j=="numTs") {
        meth@data.numTs
    } else if(j=="numCs") {
        meth@data.numCs
    } else {
        whichEle = match(j, paste("coverage", 1:NCOL(meth@data.coverage), sep=""))
        if(!is.na(whichEle)) return(meth@data.coverage[i,whichEle])
        whichEle = match(j, paste("numCs", 1:NCOL(meth@data.coverage), sep=""))
        if(!is.na(whichEle)) return(meth@data.numCs[i,whichEle])
        whichEle = match(j, paste("numTs", 1:NCOL(meth@data.coverage), sep=""))
        if(!is.na(whichEle)) return(meth@data.numTs[i,whichEle])
    }
}

#' Subset methylSigData object based on a range
#'
#' @param i Start position
#' @param j End position
#'
#' @rdname methylSigData-class
#'
#' @keywords methods internal
setMethod("[", signature(x = "methylSigData", i = "ANY", j = "ANY", drop="ANY"),
    function(x, i, j, drop) {
        if(missing(i) && missing(j)) {
            x
        } else {
            if(missing(i)) i = 1:NROW(x@data.chr)
            if(!missing(j) && is.character(j)) {
                methylSig.subDataElement(x,i,j)
            } else {
                methylSig.subData(x,i,j,drop)
            }
        }
    }
)

#' This method is for printing a methylSigData object.
#'
#' @param object A methylSigData object.
#'
#' @examples
#' show("methylSigData")
#'
#' @rdname methylSigData-class
#'
#' @keywords methods internal
setMethod("show", "methylSigData",
  function(object) {
    cat("methylSigData object with", format(NROW(object@data.start), big.mark=",",small.interval=3), "rows\n")
    cat("--------------------------\n")
    printRange = 1:min(10,NROW(object@data.coverage))
    nCol = NCOL(object@data.coverage)
    printCCT = matrix(0, nrow=NROW(printRange), ncol=3*nCol)
    printCCT[,3*(1:nCol)-2] = object@data.coverage[printRange,]
    printCCT[,3*(1:nCol)-1] = object@data.numCs[printRange,]
    printCCT[,3*(1:nCol)] = object@data.numTs[printRange,]

    colnames(printCCT) = paste(c("coverage", "numCs", "numTs"), rep(1:nCol,each=3), sep="")

    printData = data.frame(chr=object@data.chr[printRange], start=object@data.start[printRange],
                           end=object@data.end[printRange], strand=object@data.strand[printRange],
                           printCCT)

    print(printData)
    cat("--------------------------\n")
    cat("sample.ids:",object@sample.ids,"\n", sep=" ")
    cat("treatment:",object@treatment,"\n", sep=" ")
    cat("destranded:",object@destranded,"\n", sep=" ")
    cat("resolution:", object@resolution,"\n", sep=" ")
    cat("options:", object@options,"\n", sep=" ")
  }
)
