\name{methylSig.tfbsEnrichTest}
\alias{methylSig.tfbsEnrichTest}
\title{
    perform transcription factor enrichment test among differentially methylated cytosines or regions
}
\description{
    This function tests for enriched transcription binding sites among differentially methylated sites or regions using a binomial test. It also has an option to generate a bar plot to show the significance levels for the enriched transcription factor binding sites.
}
\usage{
methylSig.tfbsEnrichTest(myDiff, dmcList, tfbsInfo, plot = TRUE, 
                        max.plot.num = 10, main = "tfbs Enrichment")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myDiff}{
    methylSigDiff object that contains all CpG sites that are tested for differential methylation.
}
  \item{dmcList}{
    a vector of TRUE or FALSE with the same length as myDiff defining the DMCs or DMRs. A value of TRUE for dmcList[i] represents that the ith CpG site in myDiff is a DMC.
}
  \item{tfbsInfo}{
    an object that contains transcription factor binding sites information.
}
  \item{plot}{
    If plot = TRUE, this function draws a plot. Otherwise no plot is given. Default is TRUE.
}
  \item{max.plot.num}{
    maximum number of transcription factor binding sites on the plot. Default is 10.
}
  \item{main}{
    an overall title for the plot.
}
}
\details{
    Likelihood ratio test is used based on the binomial distribution.
}

\value{
    p-values of enriched transcription factor binding sites.
}

\seealso{
    \code{\link{getTFBSInfo}}.
}

\examples{
data(methylSig)

tfbsInfo = getTFBSInfo(system.file("annotation", "tfbsUniform.txt", 
                         package = "methylSig"))

############ tfbs Enrichment Test ###############
pvalue = methylSig.tfbsEnrichTest(myDiffSigboth,
          myDiffSigboth[,"qvalue"]<0.05 & abs(myDiffSigboth[,"meth.diff"])>25,
          tfbsInfo)  

}
\keyword{tfbsAnnotation}
