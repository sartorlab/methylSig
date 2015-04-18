\name{getCpGInfo}
\alias{getCpGInfo}
\title{
    read CpG island information from a text file. 
}

\description{
    This function reads CpG island information from a text or bed file for future annotation. 
}

\usage{
    getCpGInfo(infile)
}
\arguments{
  \item{infile}{
    a filename to read.
}
}
\details{
    This function reads CpG island informaiton for future annotation. The same genome assembly of data should be used in annotation. This function does not check genome assembly type or version.
}
\value{
     CpG island information data base that can be used in function cpgAnnotation.
}

\seealso{
    \code{\link{cpgAnnotationPlot}} and \code{\link{cpgAnnotation}}.
}

\examples{
cpgInfo = getCpGInfo(system.file("annotation", "cpgi.hg18.bed.txt", 
                       package = "methylSig"))
}

\keyword{CpGAnnotation}
