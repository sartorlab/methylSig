\name{getRefgeneInfo}
\alias{getRefgeneInfo}
\title{
    read RefGene model informaiton.
}
\description{
    This function reads RefGene model information for future annotation. 
}
\usage{
getRefgeneInfo(infile)
}

\arguments{
  \item{infile}{
    a filename that contains refGene information.
}
}

\details{
    This function reads RefGene informaiton for future annotation. Same genome assembly of data should be used in annotation. This function does not check genome assembly type nor version.
}
\value{
    refGene information object that can be used in function \code{\link{refGeneAnnotation}}.
}

\seealso{
    \code{\link{refGeneAnnotation}} and \code{\link{refGeneAnnotationPlot}}
}
\examples{
refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt", 
                             package = "methylSig"))
}
\keyword{refGeneAnnotation}
