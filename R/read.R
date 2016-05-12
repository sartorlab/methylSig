# Called by methylSigReadData
methylSigReadDataSingleFile <- function(fileIndex, fileList, header, minCount, maxCount, destranded, filterSNPs, quiet=FALSE) {
	if(!quiet) {
		message(sprintf('Reading file (%s/%s) -- %s', fileIndex, length(fileList), fileList[[fileIndex]]))
	}

    df = read.table(fileList[[fileIndex]], header=header, stringsAsFactors=FALSE)
    names(df) <- c("id", "chr", "start", "strand", "coverage", "numCs", "numTs")

	# Convert all strand information into +/-/* from {+,-,*,F,R,.}
	df$strand[which(df$strand == 'F')] = '+'
	df$strand[which(df$strand == 'R')] = '-'
	df$strand[which(df$strand == '.')] = '*'
	df$strand = factor(df$strand, levels=c('+','-','*'))

    # At this point, numCs and numTs are actually freqC and freqT from methylKit.
	# Filter sites with C + T < 95%
    freqInvalidList = which(df$numCs + df$numTs < 95)
    if(length(freqInvalidList) > 0) {
		if(!quiet) {
			message(sprintf('File (%s/%s) Sites with numCs + numTs < 95: %s/%s = %s',
				fileIndex, length(fileList), length(freqInvalidList), nrow(df), signif(length(freqInvalidList) / nrow(df), 3)))
		}
		df$coverage[freqInvalidList] = 0
    } else {
		if(!quiet) {
			message(sprintf('File (%s/%s) Sites with numCs + numTs < 95: 0 / %s = 0',
				fileIndex, length(fileList), nrow(df)))
		}
	}

	# Now numCs and numTs have frequencies replaced by counts
	df$numCs = round(df$numCs * df$coverage / 100)
	df$numTs = round(df$numTs * df$coverage / 100)

	# Destrand (or not)
	if(destranded == TRUE) {
		if(!quiet) {
			message(sprintf('File (%s/%s) Destranding',
				fileIndex, length(fileList)))
		}

		# Shift positions on - strand back by 1bp to combine with + strand
		destrandList = which(df$strand == "-" | df$strand == "R")
        df$start[destrandList] = df$start[destrandList] - 1

		# Need to re-sort the df in case destranding puts things out of order
		# For example,
		# chr1	6	+
		# chr1	6	-
		# Would be destranded to
		# chr1	6	+
		# chr1	5	+
		df = df[order(df$chr, df$start),]

		# Some ingredients for the hash / destranding
		MAXBASE = 10^{ceiling(log10(max(df$start) + 1))}
		uniqueChr = sort(unique(df$chr))

		# Hash the locations
		uniqueLoc = unique(as.numeric(as.factor(df$chr)) * MAXBASE + df$start)
		numSites = length(uniqueLoc)

		# Create empty matrices where file data will be columns
	    chrom = pos = coverage = numCs = numTs = rep.int(0, numSites)
	    strand = factor(rep(NA, numSites), levels=c('+','-','*'))

		# Determine locations of data in uniqueLoc. This should behave so that
		# base shifted - strand sites have the same "location".
		location =  findInterval( (as.numeric(as.factor(df$chr)) * MAXBASE) + df$start, uniqueLoc)

		# Deal with locations having strand NA and set to *
		# settingList = is.na(strand[location])
		# strand[location][settingList] = df$strand[settingList]
		# settingList = !settingList & (strand[location] != df$strand)
		# strand[location][settingList] = "*"
		strand = factor(rep.int('*', numSites), levels=c('+','-','*'))

		# Deal with positive strand
		forward = (df$strand == "+") | (df$strand == "F")
		chrom[location[forward]] = df$chr[forward]
		pos[location[forward]] = df$start[forward]
		coverage[location[forward]] = df$coverage[forward]
		numCs[location[forward]] = df$numCs[forward]
		numTs[location[forward]] = df$numTs[forward]

		# Deal with reverse strand
		reverse = (df$strand == "-") | (df$strand == "R")
		chrom[location[reverse]] = df$chr[reverse]
		pos[location[reverse]] = df$start[reverse]
		coverage[location[reverse]] = coverage[location[reverse]] + df$coverage[reverse]
		numCs[location[reverse]] = numCs[location[reverse]] + df$numCs[reverse]
		numTs[location[reverse]] = numTs[location[reverse]] + df$numTs[reverse]

		# Rebuild the data.frame
		df = data.frame(
			chr = chrom,
			start = pos,
			strand = strand,
			coverage = coverage,
			numCs = numCs,
			numTs = numTs,
			stringsAsFactors=F)
    } else {
		if(!quiet) {
			message(sprintf('File (%s/%s) Not destranding',
				fileIndex, length(fileList)))
		}
		# Remove the first column of the data.frame
		df = df[,2:7]
	}

	# Filter for minCount or maxCount
	countInvalidList = which(df$coverage > maxCount | df$coverage < minCount)
	if(length(countInvalidList) > 0) {
		if(!quiet) {
			message(sprintf('File (%s/%s) Sites > maxCount or < minCount: %s/%s = %s',
				fileIndex, length(fileList), length(countInvalidList), nrow(df), signif(length(countInvalidList) / nrow(df), 3)))
		}
		df$coverage[countInvalidList] = 0
	} else {
		if(!quiet) {
			message(sprintf('File (%s/%s) Sites > maxCount or < minCount: 0 / %s = 0',
				fileIndex, length(fileList), nrow(df)))
		}
	}

	# Filter C > T SNPs
	if(filterSNPs) {
		data('CT_SNPs_hg19', envir=environment())

		df_gr = GRanges(seqnames=df$chr, ranges=IRanges(start=df$start, end=df$start))
		overlaps = findOverlaps(df_gr, CT_SNPs_hg19)
		snpInvalidList = overlaps@queryHits

		if(length(snpInvalidList) > 0) {
			if(!quiet) {
				message(sprintf('File (%s/%s) Sites overlapping C > T SNP: %s/%s = %s',
					fileIndex, length(fileList), length(snpInvalidList), nrow(df), signif(length(snpInvalidList) / nrow(df), 3)))
			}
			df$coverage[snpInvalidList] = 0
		} else {
			if(!quiet) {
				message(sprintf('File (%s/%s) Sites overlapping C > T SNP: 0 / %s = 0',
					fileIndex, length(fileList), nrow(df)))
			}
		}
    } else {
		if(!quiet) {
			message(sprintf('File (%s/%s) filterSNPs == FALSE',
				fileIndex, length(fileList)))
		}
	}

	return(df[df$coverage>0,])
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

	# Collect chromosomes present in each of the files and determine the maximum
	# start location. (MAXBASE PURPOSE NOT CLEAR YET)
    MAXBASE = 0
    uniqueChr = NULL
    for(fileIndex in 1:n.files) {
         uniqueChr = c(uniqueChr, chrList[[fileIndex]]$chr)
         MAXBASE = max(MAXBASE, max(chrList[[fileIndex]]$start))
    }

	# Determine the unique chromosomes and order them
    uniqueChr = unique(uniqueChr)
    uniqueChr = uniqueChr[order(uniqueChr)]

	# Basically round MAXBASE to the nearest power of 10
    MAXBASE = 10^{ceiling(log10(MAXBASE + 1))}

	# Create a hash for unique locations
	uniqueLoc = NULL
    for(fileIndex in 1:n.files) {
        chrList[[fileIndex]]$chr = factor(chrList[[fileIndex]]$chr, levels=uniqueChr)
        uniqueLoc = unique(c(uniqueLoc, as.numeric(chrList[[fileIndex]]$chr) * MAXBASE+chrList[[fileIndex]]$start))
    }

	# Determine only the unique locations and the number to create the empty
	# matrices for aggregation
    uniqueLoc = uniqueLoc[order(uniqueLoc)]
    sizeRet = NROW(uniqueLoc)

	# Create empty matrices where file data will be columns
    coverage = numCs = numTs = matrix(0, nrow=sizeRet, ncol=n.files)
    strand = factor(rep(NA, sizeRet), levels=levels(chrList[[1]]$strand))

	# Aggregate
    for(fileIndex in 1:n.files) {
        # if(quiet == FALSE) {
        #     message(sprintf('(%s)', fileIndex))
        # }
        location =  findInterval( (as.numeric(chrList[[fileIndex]]$chr) * MAXBASE+chrList[[fileIndex]]$start) , uniqueLoc)

        strand[location] = chrList[[fileIndex]]$strand
        coverage[location,fileIndex] = chrList[[fileIndex]]$coverage
        numCs[location,fileIndex] = chrList[[fileIndex]]$numCs
        numTs[location,fileIndex] = chrList[[fileIndex]]$numTs
    }

    options = paste("maxCount=", maxCount, " & minCount=", minCount, " & filterSNPs=", filterSNPs, sep="")
    if(!is.na(assembly)) options = paste(options, " & assembly=", assembly, sep="")
    if(!is.na(context))  options = paste(options, " & context=",  context,  sep="")
    if(!is.na(pipeline)) options = paste(options, " & pipeline=", pipeline, sep="")

    methylSig.newData(data.ids=uniqueLoc, data.chr=as.factor(uniqueChr[as.integer(uniqueLoc/MAXBASE)]),
                      data.start=uniqueLoc%%MAXBASE, data.end=uniqueLoc%%MAXBASE,
                      data.strand=strand, data.coverage = coverage, data.numTs = numTs, data.numCs = numCs,
                      sample.ids=sample.ids, treatment=treatment, destranded=destranded,
                      resolution=resolution, sample.filenames=fileList,options=options)
}
