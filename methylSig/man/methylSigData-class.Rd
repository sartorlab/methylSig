\name{methylSigData-class}
\Rdversion{1.1}
\docType{class}
\alias{methylSigData-class}
\alias{[,methylSigData,ANY,ANY-method}
\alias{show,methylSigData-method}

\title{Class \code{"methylSigData"}}
\description{
    This class is for generating 'methylSigData' object.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("methylSigData", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{data.ids}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{data.chr}:}{Object of class \code{"factor"} ~~ }
    \item{\code{data.start}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{data.end}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{data.strand}:}{Object of class \code{"factor"} ~~ }
    \item{\code{data.coverage}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{data.numCs}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{data.numTs}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{sample.ids}:}{Object of class \code{"character"} ~~ }
    \item{\code{sample.filenames}:}{Object of class \code{"character"} ~~ }
    \item{\code{treatment}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{destranded}:}{Object of class \code{"logical"} ~~ }
    \item{\code{resolution}:}{Object of class \code{"character"} ~~ }
    \item{\code{options}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "methylSigData", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{show}{\code{signature(object = "methylSigData")}: ... }
	 }
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("methylSigData")
}
\keyword{classes}
\keyword{internal}
