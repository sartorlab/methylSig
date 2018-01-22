#' Plot DM status distribution in annotations
#'
#' Generates a bar plot showing the distribution of differential methylation status of loci in selected annotations.
#'
#' @param myAnnots A \code{GRanges} object resulting from running \code{methylSigAnnotation} on the result of \code{methylSigCalc}.
#' @param annotation_order A character vector which orders and subsets the annotations for the plot.
#' @param status_order A character vector which orders and subsets the DM status.
#' @param position One of \code{fill}, \code{stack}, \code{dodge}. See \code{ggplot2} documentation for details.
#' @param plot_title The plot label.
#' @param legend_title The legend label.
#' @param x_label The x-axis label.
#' @param y_label The y-axis label.
#'
#' @return A \code{ggplot2} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' # Annotate the msig_cpgs results
#' data(data, package = 'methylSig')
#'
#' # Use the genome of msig_cpgs and build annotations for CpG features
#' genome = GenomeInfoDb::genome(msig_cpgs)
#' annots = annotatr::build_annotations(genome = genome, annotations = paste(genome, c('cpgs'), sep='_'))
#'
#' # Decide what counts as differentially methylated
#' dmcList = msig_cpgs$fdr < 0.05 & abs(msig_cpgs$meth.diff) > 25
#'
#' # Annotate
#' myDiff_annotated = methylSigAnnotation(myDiff = msig_cpgs, dmcList = dmcList, annotations = annots)
#'
#' # Set the order vectors
#' cpg_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter')
#' dm_order = c('DR','DS','No DM')
#'
#' methylSigPlotStatus(myAnnots = myDiff_annotated, annotation_order = cpg_order, status_order = dm_order,
#'     position = 'fill', plot_title = 'DM Status in CpG Annots.', legend_title = 'Annotations',
#'     x_label = 'DM Status', y_label = 'Proportion')
#'
#' @export
methylSigPlotStatus = function(myAnnots, annotation_order = NULL, status_order = NULL, position = 'fill',
    plot_title = 'DM Status', legend_title = 'Annotations', x_label = 'DM Status', y_label = 'Proportion') {

    plot = annotatr::plot_categorical(
        annotated_regions = myAnnots,
        x = 'dm_status',
        fill = 'annot.type',
        x_order = status_order,
        fill_order = annotation_order,
        position = position,
        plot_title = plot_title,
        legend_title = legend_title,
        x_label = x_label,
        y_label = y_label)

    return(plot)
}

#' Data visualization function
#'
#' Generates data visualization plot of methylation data for a specified genomic interval.
#'
#'   This function offers a unique two-tiered visualization of the methylation data depending on the zoom level. For narrow regions (<1mbp) where at most 500 CpG sites have data reads, users can visualize sample-specific coverage levels and percent methylation at each site, together with group averages, significance levels and a number of genomic annotations.
#'
#' @param myAnnots A \code{GRanges} object resulting from running \code{methylSigAnnotation} on the result of \code{methylSigCalc}.
#' @param annotation_order A character vector which orders and subsets the annotations for the plot.
#' @param status_order A character vector which orders and subsets the DM status.
#' @param bin_width A vector of two numbers (from, to) to specify the region to visualize on chromosome `chr'.
#' @param plot_title The plot label.
#' @param legend_title The legend label.
#' @param x_label The x-axis label.
#' @param y_label The y-axis label.
#' @param legend_facet_label Label explaining the meaning of the gray bars in the resulting histogram.
#' @param legend_cum_label Label explaining the meaning of the red outlined bars in the resulting histogram.
#'
#' @return A \code{ggplot2} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' # Annotate the msig_cpgs results
#' data(data, package = 'methylSig')
#'
#' # Use the genome of msig_cpgs and build annotations for CpG features
#' genome = GenomeInfoDb::genome(msig_cpgs)
#' annots = annotatr::build_annotations(genome = genome, annotations = paste(genome, c('cpgs'), sep='_'))
#'
#' # Decide what counts as differentially methylated
#' dmcList = msig_cpgs$fdr < 0.05 & abs(msig_cpgs$meth.diff) > 25
#'
#' # Annotate
#' myDiff_annotated = methylSigAnnotation(myDiff = msig_cpgs, dmcList = dmcList, annotations = annots)
#'
#' # Set the order vectors
#' cpg_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter')
#' status_order = c('hyper','hypo','none')
#'
#' methylSigPlotDiff(myAnnots = myDiff_annotated, annotation_order = cpg_order, status_order = status_order,
#'     bin_width = 10, plot_title = 'Meth. Diff. in CpG Annots.', x_label = 'DM Status', y_label = 'Proportion',
#'     legend_facet_label = 'Meth. Diff. in Annotation', legend_cum_label = 'Meth. Diff. Overall')
#'
#' @export
methylSigPlotDiff = function(myAnnots, annotation_order = NULL, status_order = NULL, bin_width = 10,
    plot_title = 'Methylation Differences', x_label = 'Methylation Difference', y_label = 'Density',
    legend_facet_label = 'Meth. Diff. in Annotation', legend_cum_label = 'Meth. Diff. Overall') {

    facet_order = list(annotation_order, status_order)

    plot = annotatr::plot_numerical(
        annotated_regions = myAnnots,
        x = 'meth.diff',
        facet = c('annot.type', 'dm_status'),
        facet_order = facet_order,
        bin_width = bin_width,
        plot_title = plot_title,
        x_label = x_label,
        y_label = y_label,
        legend_facet_label = legend_facet_label,
        legend_cum_label = legend_cum_label)

    return(plot)
}

#' Data visualization function
#'
#' Generates data visualization plot of methylation data for a specified genomic interval.
#'
#'   This function offers a unique two-tiered visualization of the methylation data depending on the zoom level. For narrow regions (<1mbp) where at most 500 CpG sites have data reads, users can visualize sample-specific coverage levels and percent methylation at each site, together with group averages, significance levels and a number of genomic annotations.
#'
#' @param meth A \code{link{methylSigData-class}} object.
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
