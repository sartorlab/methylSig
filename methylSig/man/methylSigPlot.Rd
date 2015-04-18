\name{methylSigPlot}
\alias{methylSigPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   data visualization function.
}
\description{
   This function generates data visualization plot of methylation data for a specified genomic interval. 
}
\usage{
methylSigPlot(meth,chr,loc.range,groups=c(Treatment=1,Control=0), 
              cpgInfo, refGeneInfo, myDiff, tfbsInfo, 
              noGroupEst = FALSE, noDataLine = TRUE, cex = 1, 
              tfbsDense = TRUE, sigQ = 0.05)
}
\arguments{
  \item{meth}{
    a `methylSigData' object.
}
  \item{chr}{
    chromosome in character that matches in methylation data (meth).
}
  \item{loc.range}{
    a vector of two numbers to specify the region to visualize on chromosome `chr'. 
}
  \item{groups}{
    a vector of two numbers to specify two groups to show in the plot.
}
  \item{cpgInfo}{
    CpG island information to annotate. If missing, no CpG island information will be shown.
}
  \item{refGeneInfo}{
    refGene information to annotate. If missing, no refGene island information will be shown.
}
  \item{myDiff}{
    differential methylation analysis results from `methylSigCalc'. If missing, no ``-log10[pvalues]'' will be shown.
}
  \item{tfbsInfo}{
    transcription factor binding site informaiton from `getTFBSInfo' function. If missing, no tfbs will be shown.
}
  \item{noGroupEst}{
    if noGroupEst = TRUE, group estimates will not be drawn in the plot. Default is FALSE
}
  \item{noDataLine}{
    if noDataLine = FALSE, a line will link the same letter that represents the same sample. Default is TRUE.
}
  \item{cex}{
    A numerical value giving the amount by which plotting letters of data should be magnified relative to the default.
}
  \item{tfbsDense}{
    if tfbsDense = TRUE, the plot of transcription factor binding sites is using dense mode. If the region is broader, the tfbsDense is automatically set to TRUE.
}
  \item{sigQ}{
    maximum significant q-value to be drawn as red line. Default is 0.05.
}
}
\details{
    This function offers a unique two-tiered visualization of the methylation data depending on the zoom level. For narrow regions (<1mbp) where at most 500 CpG sites have data reads, users can visualize sample-specific coverage levels and percent methylation at each site, together with group averages, significance levels and a number of genomic annotations.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(methylSig)

cpgInfo = getCpGInfo(system.file("annotation", "cpgi.hg18.bed.txt",
                  package = "methylSig"))
tfbsInfo = getTFBSInfo(system.file("annotation", "tfbsUniform.txt",
                  package = "methylSig"))
refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt",
                  package = "methylSig"))

############ broader region ####################
methylSigPlot(meth, "chr21", c(43000000, 44000000), groups=c(1,0), 
              cpgInfo=cpgInfo, refGeneInfo=refGeneInfo, 
              myDiff=myDiffSigboth, tfbsInfo=tfbsInfo,
              tfbsDense=FALSE, sigQ=0.05)

############ narrower region ####################
methylSigPlot(meth, "chr21", c(43800000, 43900000), groups=c(1,0), 
              cpgInfo=cpgInfo, refGeneInfo=refGeneInfo, 
              myDiff=myDiffSigboth, tfbsInfo=tfbsInfo,
              tfbsDense=FALSE, sigQ=0.05)
}
