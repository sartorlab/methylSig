\name{refGeneAnnotation}
\alias{refGeneAnnotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    annotation function using RefGene models
}
\description{
    This function annotates the CpG sites to promoter, coding DNA sequence (CDS), exons, 3' untranslated region, 5' untranslated region, nocoding RNA region, or integenic regions using the RefGene gene models.
}
\usage{
refGeneAnnotation(refGeneInfo, myDiff)
}

\arguments{
  \item{refGeneInfo}{
     refGene information object from function \code{\link{getRefgeneInfo}}.
}
  \item{myDiff}{
     `methylSigDiff' object to be used for annotation.
}
}

\details{
     This function generates an object to annotate refGene information. The same genome assembly should be used for `myDiff' and `refGeneInfo'.

     The promoter region is defined as 1,000bp upstream from a transcription starting site(TSS) to the TSS.
}
\value{
     A logical matrix of 6 rows and number of columns equal to the size of `myDiff' with each entry indicating whether the related sites in `myDiff' are in a promoter, CDS, 3' UTR, 5' UTR, noncoding RNA or exon regions. Regions can overlap because of a single gene can have multiple isoforms.
}

\seealso{
    \code{\link{refGeneAnnotationPlot}} and \code{\link{getRefgeneInfo}}
}

\examples{
data(methylSig)
refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt", 
                                package = "methylSig"))

refGeneAnn = refGeneAnnotation(refGeneInfo,myDiffSigboth)
}
    
\keyword{ refGeneAnnotation }
