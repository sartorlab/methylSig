\name{methylSigDiff-class}
\Rdversion{1.1}
\docType{class}
\alias{methylSigDiff-class}
\alias{[,methylSigDiff,ANY,missing-method}
\alias{[,methylSigDiff,ANY,character-method}
\alias{show,methylSigDiff-method}

\title{Class \code{"methylSigDiff"}}
\description{
    This class is to generate 'methylSigDiff' object.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("methylSigDiff", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{data.ids}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{data.chr}:}{Object of class \code{"factor"} ~~ }
    \item{\code{data.start}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{data.end}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{data.strand}:}{Object of class \code{"factor"} ~~ }
    \item{\code{results}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{sample.ids}:}{Object of class \code{"character"} ~~ }
    \item{\code{treatment}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sample.filenames}:}{Object of class \code{"character"} ~~ }
    \item{\code{destranded}:}{Object of class \code{"logical"} ~~ }
    \item{\code{resolution}:}{Object of class \code{"character"} ~~ }
    \item{\code{options}:}{Object of class \code{"character"} ~~ }
    \item{\code{data.options}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "methylSigDiff", i = "ANY", j = "character")}: ... }
    \item{[}{\code{signature(x = "methylSigDiff", i = "ANY", j = "missing")}: ... }
    \item{show}{\code{signature(object = "methylSigDiff")}: ... }
	 }
}

\examples{
showClass("methylSigDiff")
}
\keyword{internal}
\keyword{classes}
