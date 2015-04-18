\name{methylSigTile}
\alias{methylSigTile}
\title{
    obtain tiled methylation data in non-overlapping continuous windows
}
\description{
    This funciton tiles data within windows of a given width across genome.  It gives total number of methylated cytosines and total number of reads (coverage) within each tiled window.  This is used prior to the \code{\link{methylSigCalc}} function when the user prefers to conduct a tiled analysis instead of a base specific analysis. Tiling may provide higher power to detect significant differences, especially for experiments with low coverage.
}
\usage{
methylSigTile(meth, win.size = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meth}{
    a `methylSigData' object used to tile data.
}
  \item{win.size}{
    an integer value indicating the desired window size in bps. Default is 25.
}
}
\value{
   `methylSigData' object. 
}

\examples{
data(methylSig)

methTile = methylSigTile(meth)
}
