\name{internalFunctions}
\alias{methylSig_weightFunc}
\alias{methlSig.subDataElement}
\alias{methlSig.subDiffElement}
\alias{methylSig.newData}
\alias{methylSig.newDiff}
\alias{cpgAnnotationByChrom}
\alias{getTFBSCountByChrom}
\alias{is.TFBSByChrom}
\alias{refGeneAnnotationByChrom}
\alias{methylSig_derivativeMu}
\alias{methylSig_derivativePhi}
\alias{methylSig_dataProcess}
\alias{methylSigReadDataSingleFile}
\alias{methylSig_logLik}
\alias{methylSig.subData}
\alias{methylSig.subDiff}
\alias{methylSigDf}

\title{
    default weight funciton.
    function called from subsetting methylSigData.
    function called from subsetting methylSigDiff.
    creating new methylSigData.
    creating new methylSigDiff.
    get CpG annotation information by chromosome
    get tfbs annotation results within a chromosome
    test whether the CpG sites are within transcription factor binding peaks.
    Annotation function using RefGene model for the sites in the same chromosome
}
\description{
    internal functions.
}
\usage{
methylSig_weightFunc(u)
methlSig.subDataElement(meth, i, j)
methlSig.subDiffElement(meth, i, j)
methylSig.newData(data.ids, data.chr, data.start, data.end, data.strand, data.coverage, data.numCs, data.numTs, sample.ids = "", sample.filenames = "", treatment, destranded = TRUE, resolution = "", options = "")
methylSig.newDiff(data.ids, data.chr, data.start, data.end, data.strand, results, sample.ids = "", sample.filenames = "", treatment, destranded = TRUE, resolution = "", options = "", data.options = "")
}

\arguments{
  \item{u}{
    input value range from -1 to 1.
}
}

\keyword{internal}

