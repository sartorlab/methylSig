# Called by methylSigReadData
methylSigReadDataSingleFile <- function(fileIndex, fileList, header, minCount, maxCount, destranded, filterSNPs, quiet=FALSE) {
    if(quiet==FALSE) {
      message(sprintf('Reading file (%s / %s) -- %s', fileIndex, nrow(fileList), fileList[[fileIndex]]))
    }
    chr = read.table(fileList[[fileIndex]], header=header, stringsAsFactors=FALSE)
    #### order for base ####
    ##chr = chr[order(chr$base),]
    ####
    names(chr) <- c("id",  "chr", "start", "strand", "coverage", "numCs", "numTs")

    # At this point, numCs and numTs are actually freqC and freqT from methylKit.
    freqInvalidList = which(chr$numCs + chr$numTs < 95)
    if(length(freqInvalidList) > 0) {
        message(sprintf('File (%s / %s) Sites with numCs + numTs < 95: %s / %s = %s',
          fileIndex, nrow(fileList), nrow(freqInvalidList), nrow(chr), signif(nrow(freqInvalidList) / nrow(chr), 3)))
        chr$coverage[freqInvalidList] = 0
    }

    # Now numCs and numTs have frequencies replaced by counts
    chr$numCs<- round(chr$numCs * chr$coverage / 100)
    chr$numTs<- round(chr$numTs * chr$coverage / 100)

    # countInvalidList = which(chr$coverage > maxCount | chr$coverage < minCount)
    # message(sprintf('File (%s / %s) Sites > maxCount or < minCount: %s / %s = %s',
    #   fileIndex, nrow(fileList), nrow(countInvalidList), nrow(chr), signif(nrow(countInvalidList) / nrow(chr), 3)))
    # chr$coverage[countInvalidList] = 0

    if(destranded == TRUE) {
        destrandList = which(chr$strand == "-" | chr$strand == "R")
        chr$start[destrandList] = chr$start[destrandList] - 1
    }

    if(filterSNPs) {
        data('CT_SNPs_hg19', envir=environment())
        chr_gr = GRanges(seqnames=chr$chr, ranges=IRanges(start=chr$start, end=chr$start))

        overlaps = findOverlaps(chr_gr, CT_SNPs_hg19)
        snpInvalidList = overlaps@queryHits

        message(sprintf('File (%s / %s) Sites overlapping C > T SNP: %s / %s = %s',
          fileIndex, nrow(fileList), nrow(snpInvalidList), nrow(chr), signif(nrow(snpInvalidList) / nrow(chr), 3)))

        chr$coverage[snpInvalidList] = 0
    }

    #chr$chr = as.factor(chr$chr)
    chr$strand = as.factor(chr$strand)
    ## When only F observed in strand column, as.factor convert F to FALSE
    levels(chr$strand) = list("+"="F","-"="R","*"="*", "+"="FALSE")

    chr[chr$coverage>0,2:7]
}

#' Read methylation score files to make a 'methylSigData' object.
#'
#' This function reads methylation score files (having columns chrBase, chr, base, strand, coverage, freqC, freqT) to make a \code{methylSigData-class} object that can be used in differential methylation analysis.
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
#' @return A \code{methylSigData-class} object.
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
        for(i in 1:n.files) {
          chrList[[i]] =  methylSigReadDataSingleFile(i, fileList, header=header, minCount=minCount, maxCount=maxCount, destranded, filterSNPs, quiet)
        }
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
            message(sprintf('(%s)', fileIndex))
        }
        location =  findInterval( (as.numeric(chrList[[fileIndex]]$chr) * MAXBASE+chrList[[fileIndex]]$start) , uniqueLoc)
        if(destranded == FALSE) {
            strand[location] = chrList[[fileIndex]]$strand
            coverage[location,fileIndex] = chrList[[fileIndex]]$coverage
            numCs[location,fileIndex] = chrList[[fileIndex]]$numCs
            numTs[location,fileIndex] = chrList[[fileIndex]]$numTs
        } else {
            # Deal with locations having strand NA and set to *
            settingList = is.na(strand[location])
            strand[location][settingList] = chrList[[fileIndex]]$strand[settingList]
            settingList = !settingList & (strand[location] != chrList[[fileIndex]]$strand)
            strand[location][settingList] = "*"

            # Deal with positive strand
            forward = (chrList[[fileIndex]]$strand == "+")
            coverage[location[forward],fileIndex] = chrList[[fileIndex]]$coverage[forward]
            numCs[location[forward],fileIndex] = chrList[[fileIndex]]$numCs[forward]
            numTs[location[forward],fileIndex] = chrList[[fileIndex]]$numTs[forward]

            # Deal with reverse strand
            reverse = (chrList[[fileIndex]]$strand == "-")
            coverage[location[reverse],fileIndex] = coverage[location[reverse],fileIndex] + chrList[[fileIndex]]$coverage[reverse]
            numCs[location[reverse],fileIndex] = numCs[location[reverse],fileIndex] + chrList[[fileIndex]]$numCs[reverse]
            numTs[location[reverse],fileIndex] = numTs[location[reverse],fileIndex] + chrList[[fileIndex]]$numTs[reverse]
        }
    }

    # Filter by minCount and maxCount while also modifying chr, uniqueLoc, strand, numCs, numTs
    # This is in the idiom of the what's been implemented
      countInvalidList = NULL
      for(fileIndex in 1:n.files) {
        # Record the countInvalidList
        countInvalidList = unique( c(countInvalidList, which(
          ((coverage[, fileIndex] != 0) & (coverage[, fileIndex] < minCount)) |
          ((coverage[, fileIndex] != 0) & (coverage[, fileIndex] > maxCount)) ) ) )
      }

      uniqueChr = uniqueChr[-countInvalidList]
      uniqueLoc = uniqueLoc[-countInvalidList]
      strand = strand[-countInvalidList]
      coverage = coverage[-countInvalidList, ]
      numCs = numCs[-countInvalidList, ]
      numTs = numTs[-countInvalidList, ]

      message(sprintf('Sites > maxCount or < minCount: %s / %s = %s',
        length(countInvalidList), nrow(coverage) + length(countInvalidList),
        signif( length(countInvalidList) / (nrow(coverage) + length(countInvalidList)) , 3)))

#    strand = as.factor(strand)
#    levels(strand) = list("+"="1","-"="2","*"="3")

    options = paste("maxCount=", maxCount, " & minCount=", minCount, " & filterSNPs=", filterSNPs, sep="")
    if(!is.na(assembly)) options = paste(options, " & assembly=", assembly, sep="")
    if(!is.na(context))  options = paste(options, " & context=",  context,  sep="")
    if(!is.na(pipeline)) options = paste(options, " & pipeline=", pipeline, sep="")

    # numTs[coverage==0] = NA
    # numCs[coverage==0] = NA
    # coverage[coverage==0] = NA
    methylSig.newData(data.ids=uniqueLoc, data.chr=as.factor(uniqueChr[as.integer(uniqueLoc/MAXBASE)]),
                      data.start=uniqueLoc%%MAXBASE, data.end=uniqueLoc%%MAXBASE,
                      data.strand=strand, data.coverage = coverage, data.numTs = numTs, data.numCs = numCs,
                      sample.ids=sample.ids, treatment=treatment, destranded=destranded,
                      resolution=resolution, sample.filenames=fileList,options=options)
}
