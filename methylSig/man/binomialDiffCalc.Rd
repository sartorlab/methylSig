\name{binomialDiffCalc}
\alias{binomialDiffCalc}
\title{
    differential methylation analysis using binomial model
}

\description{
    This function calculates differential methylation statistics using a binomial-based approach. See `Warning' message below.
}
\usage{
binomialDiffCalc(meth, groups = c(Treatment = 1, Control = 0), 
                  min.per.group = c(3, 3))
}
\arguments{
  \item{meth}{
    a `methylSigData' object to calculate differential methylation statistics. It can be obtained using `methylSigReadData'.
}
  \item{groups}{
    a vector of two numbers specify two groups to compare. See `treatment' argument of \code{\link{methylSigReadData}} function. Default is c(Treatment=1,Control=0). 
}
  \item{min.per.group}{
    a vector with two numbers that specify the minimum numbers of samples required to be qualify as defferentially methylation region.  If it is a single number, both groups will use it as the minimum requried number of samples. Default is c(3,3).
}
}
\details{
    This function uses a binomial-based model to calculate differential methylation statistics. It is nearly identical to the `calculateDiffMeth' function in the methylKit R package except that only likelihood ratio test and p.adjust() with method=``BH'' are used to calculate significance levels. It is significantly faster than `calculateDiffMeth' function.
}

\value{
  `methylSigDiff' object that contains the differential methylation statistics and chromosomal locations. p.adjust with method="BH" option is used for P-value correction.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{
    Warning: This function does not take into account the variability among samples in each group being compared.
} 

\seealso{
    \code{\link{methylSigCalc}}
}
\examples{
data(methylSig)

myDiff = binomialDiffCalc(meth)
}
\keyword{differentialMethylation}
