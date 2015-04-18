\name{cpgAnnotation}
\alias{cpgAnnotation}
\title{
CpG island annotation function
}
\description{
This function annotates the CpG sites or regions to CpG island, CpG shore, CpG shelf or interCGI regions.
}
\usage{
cpgAnnotation(cpgInfo, myDiff)
}

\arguments{
  \item{cpgInfo}{
     a CpG island information object.
}
  \item{myDiff}{
     a `methylSigDiff' object to be used for annotation.
}
}

\details{
This function generates an object to annotate CpG island information. The same genome assembly should be used for `myDiff' and `cpgInfo'. 

CpG shores are defined as the region outside CpG islands but within 2,000 bp of any CpG island. CpG shelves are defined as the region within 2,000 bp away from a CpG shore.  When regions overlap, the priority is CpG island, then CpG shore. Regions >4,000 bp from a CpG island are deemed interCGI.
}
\value{
A factor vector with the same size as `myDiff' that indicates whether the related sites in `myDiff' are in C pG island, CpG shore, CpG shelf or interCGI region. 
}

\seealso{
\code{\link{cpgAnnotationPlot}} and \code{\link{getCpGInfo}}.
}
\examples{
data(methylSig)

cpgInfo = getCpGInfo(system.file("annotation", "cpgi.hg18.bed.txt",
                                  package = "methylSig"))

cpgAnn = cpgAnnotation(cpgInfo,myDiffSigboth)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{CpGAnnotation}
