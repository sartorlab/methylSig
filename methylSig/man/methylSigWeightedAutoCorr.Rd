\name{methylSigWeightedAutoCorr}
\alias{methylSigWeightedAutoCorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    weighted auto correlation of methylation rates
}
\description{
    This funciton calculates the weighted auto correlation of methylation rates based on the coverage level at each CpG site.
}
\usage{
methylSigWeightedAutoCorr(chr, start, methRates, 
          coverage, lags = 2)
}
\arguments{
  \item{chr}{
    a single string indicating chromosome from `methylSigData' or `methylSigDiff' objects.
}
  \item{start}{
    a vector of CpG site loci from `methylSigData' or `methylSigDiff' objects.
}
  \item{methRates}{
    a numeric vector of methylation rates at each CpG site.
}
  \item{coverage}{
    a numeric vector of coverage levels (number of reads) at each CpG site.
}
  \item{lags}{
    autocorrelation lags. This argument can be a vector if you want to calculate multiiple lags. Default is 2.
}
}
\details{
    The weight for the locus i and j is coverage[i]*coverage[j]/(coverage[i]+coverage[j]). 
}
\value{
    auto correlation value or vector for the lags entered.
}

\examples{
data(methylSig)

lags = c(2,5,10,50)

autoCorr = methylSigWeightedAutoCorr(meth[,"chr"], meth[,"start"], 
             meth[,"numCs1"]/meth[,"coverage1"], meth[,"coverage1"],
             lags=lags)

autoCorr
}
