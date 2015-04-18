\name{methylSigCalc}
\alias{methylSigCalc}
\title{
calculates differential methylation statistics using a Beta-binomial approach
}
\description{
The function calculates differential methylation statistics between two groups of samples. This is the main function of the methylSig package, and the method most users should use to test for DMCs or DMRs. The function uses a Beta-binomial approach to calculate differential methylation statistics, accounting for variation among samples within each group.
}
\usage{
methylSigCalc(meth, groups = c(Treatment = 1, Control = 0),
                 dispersion = "both", local.disp = FALSE, winsize.disp
                 = 200, local.meth = FALSE, winsize.meth = 200,
                 min.per.group = c(3, 3), weightFunc =
                 methylSig_weightFunc, T.approx = TRUE, num.cores = 1)
}
\arguments{
  \item{meth}{
    a `methylSigData' object to calculate differential methylation statistics. It can be obtained using `methylSigReadData'.
}
  \item{groups}{
    a vector of two numbers specify two groups to compare. See `treatment' argument of \code{\link{methylSigReadData}} function. Default is c(Treatment=1,Control=0).
}
  \item{dispersion}{
    a value indicating which group or groups are used to estimate the variance (dispersion). If groups are defined as c(Treatment=1,Control=0), dispersion can take values "Treatment", "Control", 1, 0 or "both". Default is "both".
}
  \item{local.disp}{
    a logical value indicating whether to use local information to improve dispersion parameter estimation. Default is False.
}
  \item{winsize.disp}{
    a number to specify the window size in basepairs for local dispersion estimation. The dispersion (variance) parameter of the groups at the particular location LOC is calculated using information from LOC - winsize.disp to LOC + winsize.disp. This argument is only activated when local.disp=TRUE. Default is 200. 
}
  \item{local.meth}{
    a logical value indicating whether to use local information to improve methylation level estimation. Default is False.
}
  \item{winsize.meth}{
    a number to specify the window size in basepairs for local methylation level estimation. The group methylation level at the particular location LOC is calculated using information from LOC - winsize.meth to LOC + winsize.meth. This argument is only activated when local.meth=TRUE. Default is 200. 
}
  \item{min.per.group}{
    a vector with two numbers that specify the minimum numbers of samples required to be qualify as defferentially methylated region.  If it is a single number, both groups will use it as the minimum requried number of samples. Default is c(3,3).
}
  \item{weightFunc}{
    a weight kernel function. The input of this function is from -1 to 1. The default is the tri-weight kernel function defined as function(u) = (1-u^2)^3. Function value and range of parameter for weight function should be from 0 to 1. 
}
  \item{T.approx}{
    a logical value indicating whether to use squared t approximation for the likelihood ratio statistics. Chi-square approximation (T.approx = FALSE) is recommended when the sample size is large.  Default is TRUE.
}
  \item{num.cores}{
    an integer denoting how many cores should be used for differential methylation calculations (only can be used in
          machines with multiple cores).
}
}
\details{
The function calculates differential methylation statistics between two groups of samples. The function uses Beta-binomial approach to calculate differential methylation statistics, accounting for variation among samples within each group. Users who wish to tile their data and test for differentially methylated regions (DMRs) instead DMCs should first use the \code{\link{methylSigTile}} function before using this function.
}
\value{
    `methylSigDiff' object containing the differential methylation statistics and locations. p.adjust with method="BH" option is used for P-value correction.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{methylSigPlot}}, \code{\link{methylSigReadData}}
}
\examples{
data(methylSig)

myDiffSig = methylSigCalc (meth)

### calculate differential methylation statistics using
### treatment group 0 to evaluate dispersion parameter.
### Also use local information to improve variance estimation.

myDiffSig = methylSigCalc (meth, dispersion=0, local.disp=TRUE, 
                           min.per.group=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ differentialMethylation }
