\name{write.methylSigDiff}
\alias{write.methylSigDiff}
\title{
    write methylSigDiff object as a text file
}
\description{
    This funciton write methylSigDiff as a text file. All options for write.table are avaiable for this function.
}
\usage{
write.methylSigDiff(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a `methylSigDiff' object.
}
  \item{...}{
    arguments that are used in function write.table().
}
}

\examples{
data(methylSig)

write.methylSigDiff(myDiffSigboth, file="myResult.txt", row.names=FALSE, quote=FALSE, sep="\t")
}
