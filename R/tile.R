#' Obtain tiled methylation data in non-overlapping continuous windows.
#'
#' This funciton tiles data within windows of a given width across genome. It gives total number of methylated cytosines and total number of reads (coverage) within each tiled window. This is used prior to the \code{\link{methylSigCalc}} function when the user prefers to conduct a tiled analysis instead of a base specific analysis for differential methylation. Tiling may provide higher power to detect significant differences, especially for experiments with low coverage.
#'
#' @param meth A \code{methylSigData-class} object used to tile data.
#' @param tiles A data.frame including columns with chr, start, end columns for predefined tiles. Those CpG sites not belonging any tile will be removed from tiled data. If \code{tiles} is not \code{NULL}, then \code{win.size} is ignored.
#' @param win.size An integer value indicating the desired window size in bps. Default is 25. If \code{tiles} is not \code{NULL}, then \code{win.size} is ignored.
#'
#' @return A \code{methylSigData-class} object.
#'
#' @examples
#' data(sampleData)
#' methTile = methylSigTile(meth)
#' @export
methylSigTile <- function(meth, tiles = NULL, win.size=25) {
   if(meth@resolution == "region") stop("Object has already been tiled")

   if(is.null(tiles)) {
     message(sprintf('Tiling by %s bp windows', win.size))
       MAXBASE = max(meth@data.start) + win.size + 1
#        MAXBASE10 = MAXBASE = 10^{ceiling(log10(MAXBASE + 1))}

       startList = as.numeric(meth@data.chr)*MAXBASE + meth@data.start
       uniqueStartList = unique(startList - (startList - 1) %% win.size)

       whichRow = match(startList - (startList-1) %% win.size, uniqueStartList)
       coverage = numTs = numCs = matrix(0,nrow=NROW(uniqueStartList), ncol=NCOL(meth@data.coverage))

       for(i in 1:NROW(startList)) {
           coverage[whichRow[i],] = coverage[whichRow[i],] + meth@data.coverage[i,]
           numCs[whichRow[i],] = numCs[whichRow[i],] + meth@data.numCs[i,]
           numTs[whichRow[i],] = numTs[whichRow[i],] + meth@data.numTs[i,]
       }

       # Need to deal with the very last end
       data.MAXBASE = uniqueStartList%%MAXBASE + win.size - 1
       data.MAXBASE[length(data.MAXBASE)] = data.MAXBASE[length(data.MAXBASE)] - win.size + 2

      methylSig.newData(data.ids=uniqueStartList, data.chr=as.factor(levels(meth@data.chr)[as.integer(uniqueStartList/MAXBASE)]),
                        data.start=uniqueStartList%%MAXBASE, data.end=data.MAXBASE,
                        data.strand=factor(rep("*",NROW(uniqueStartList))), data.coverage = coverage, data.numTs = numTs, data.numCs = numCs,
                        sample.ids=meth@sample.ids, treatment=meth@treatment, destranded=meth@destranded,
                        resolution="region", sample.filenames=meth@sample.filenames,options=meth@options)
   } else {
     message('Tiling by regions')

     # Need to have matrices corresponding to the number of samples
     numCs = numTs = matrix(0, nrow = nrow(tiles), ncol = ncol(meth@data.numCs))
     for(chr in levels(meth@data.chr)) {
         whichInTiles = which(tiles$chr == chr)
         if(length(whichInTiles) > 0) {
             # Determine the start and end locations for the tiles on this chromosome
             startList = tiles$start[whichInTiles]
             endList = tiles$end[whichInTiles]

             # Extract rows of meth are relevant
             whichVlist = which(meth@data.chr==chr)

             # Find which
             whichSTART = findInterval(startList, meth@data.start[whichVlist]+1)
             whichEND = findInterval(endList, meth@data.start[whichVlist])

             for(i in which(whichEND > whichSTART)) {
                 numCs[whichInTiles[i],] = numCs[whichInTiles[i], ] + colSums(meth@data.numCs[ whichVlist[whichSTART[i]:whichEND[i]], ], na.rm=T)
                 numTs[whichInTiles[i],] = numTs[whichInTiles[i], ] + colSums(meth@data.numTs[ whichVlist[whichSTART[i]:whichEND[i]], ], na.rm=T)
             }
         }
     }
     methylSig.newData(data.ids=1:NROW(tiles), data.chr=factor(tiles$chr), data.start=tiles$start, data.end=tiles$end,
                       data.strand=factor(rep("*",nrow(tiles))), data.coverage = numCs + numTs, data.numTs = numTs, data.numCs = numCs,
                       sample.ids=meth@sample.ids, treatment=meth@treatment, destranded=meth@destranded,
                       resolution="region", sample.filenames=meth@sample.filenames,options=meth@options)
   }
}

#' Obtain tiled methylation data by tiling (pooling) data that a particular transcription factor (TF) is predicted to bind.
#'
#' Gives total number of methylated cytosines and total number of reads (coverage) for each TF. This is used prior to the \code{\link{methylSigCalc}} function when the user prefers to conduct tests to identify significant level of hypermethylation or hypomethylation across their binding sites.
#'
#' @param meth A \code{methylSigData-class} object used to tile data.
#' @param tfbsInfo An object that contains transcription factor binding sites information.
#'
#' @return A \code{methylSigData-class} object.
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
