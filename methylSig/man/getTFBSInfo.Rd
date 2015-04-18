\name{getTFBSInfo}
\alias{getTFBSInfo}
\title{
    read transcription factor binding site information
}
\description{
    This function reads the transcription factor binding site information for future annotation. 
}
\usage{
getTFBSInfo(infile)
}
\arguments{
  \item{infile}{
    a filename to read.
}
}

\details{
    This function reads transcription factor binding sites information for future annotation. The same genome assembly of data should be used in annotation. This function does not check genome assembly type or version.
}
\value{
    an object containing transcription factor binding sites information is returned. It can be used in function \code{\link{methylSig.tfbsEnrichTest}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{methylSig.tfbsEnrichTest}},
}
\examples{
    tfbsInfo = getTFBSInfo(system.file("annotation", "tfbsUniform.txt", 
                        package = "methylSig"))
}
\keyword{tfbsAnnotation}
