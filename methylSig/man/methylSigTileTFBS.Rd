\name{methylSigTileTFBS}
\alias{methylSigTileTFBS}
\title{
    obtain tiled methylation data by tiling (pooling) data that a particular transcription factor (TF) is predicted to bind.
}
\description{
    It gives total number of methylated cytosines and total number of reads (coverage) for each TF. This is used prior to the \code{\link{methylSigCalc}} function when the user prefers to conduct tests to identify significant level of hypermethylation or hypomethylation across their binding sites.
}
\usage{
methylSigTileTFBS(meth, tfbsInfo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meth}{
    a `methylSigData' object used to tile data.
}
  \item{tfbsInfo}{
    an object that contains transcription factor binding sites information.
}
}
\value{
   `methylSigData' object. 
}
\seealso{
    \code{\link{getTFBSInfo}} and \code{\link{methylSigTile}}.
}

\examples{
data(methylSig)

methTile = methylSigTile(meth)
}
