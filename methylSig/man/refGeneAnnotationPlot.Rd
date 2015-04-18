\name{refGeneAnnotationPlot}
\alias{refGeneAnnotationPlot}

\title{
    refGene annotation plot
}
\description{
    This function generates a refGene annotation bar chart.
}
\usage{
refGeneAnnotationPlot(listFrom, main = "Plot", 
     priority = c("cds", "promoter", "noncoding", "5'utr", "3'utr"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listFrom}{
    object from fuction `refGeneAnnotation'.
}
  \item{main}{
    an overall title for the plot. Default is "Plot".
}
  \item{priority}{
    character vector indicating the priority order for assigning the annotated CpG site when overlapping occurs.  Default is coding sequence exon, promoter, noncoding exon, 5' UTR exon, 3' UTR exon, intron and then integenic region.  Only "cds", "promoter", "noncoding", "5'utr" and "3'utr" are acceptable.
}
}

\details{
    Promoter is defined as the region 1,000bp upstream from a transcription start site (TSS) to the TSS. When regions overlap, the priority to assign regions for the CpG sites is defined by `priority' argument. 
}

\seealso{
   \code{\link{refGeneAnnotation}} and \code{\link{getRefgeneInfo}}.
}

\examples{
data(methylSig)

refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt",
                  package = "methylSig"))

refGeneAnn = refGeneAnnotation(refGeneInfo, myDiffSigboth)

refGeneAnnotationPlot(refGeneAnn,main="ALL", priority=c("promoter","cds",
                            "noncoding", "5'utr", "3'utr"))
}
