\name{methylSigReadData}
\alias{methylSigReadData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   read methylation score files to make a 'methylSigData' object.
}
\description{
   This function reads methylation score files to make a 'methySigData' object that can be used in differential methylation analysis.
}

\usage{
methylSigReadData(fileList, sample.ids, assembly=NA, pipeline=NA, 
                  header=TRUE, context=NA, resolution = "base", 
                  treatment, destranded=TRUE, maxCount = 500,
                  minCount = 10, num.cores = 1, quiet = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileList}{
      vector of files to be read. 
}
  \item{sample.ids}{
      vector of sample ids.
}
  \item{assembly}{
      character string indicating the genome assembly, such as "hg18", "hg19", "mm9", or "mm10".
}
  \item{pipeline}{
      character string indicating the pipepline name that generated the data, for example, "bismark".
}
  \item{header}{
     a logical value indicating whether the score files have header or not.  Default is True.
}
  \item{context}{
    methylation context string such "CpG","CpH",or "CHH".
}
  \item{resolution}{
     a string indicating whether the input data are  base-pair or regional resolution. Either "base" or "region" is allowed. Default is "base"}
  \item{treatment}{
    a numeric vector contatining numbers to distinguish the group classification.}
  \item{destranded}{
    a logical value indicating whether destrand the reverse to forward strand. If TRUE, the reads from both will be combined. Default is TRUE. 
}
  \item{maxCount}{
    a number indicating the maximum coverage count to be included.
}
  \item{minCount}{
    a number indicating the minimum coverage count to be included.
}
  \item{num.cores}{
    number of cores to be used in reading files. Default is 1.
}
  \item{quiet}{
    a logical value. If quiet=TRUE, then this function does not show progress information. Default = FALSE.
}
}
\value{
    'methylSigData' object.
}

\seealso{
    \code{\link{methylSigCalc}}.
}
\examples{

fileList = c( system.file("extdata", "AML_1.txt", package = "methylSig"),
              system.file("extdata", "AML_2.txt", package = "methylSig"),
              system.file("extdata", "AML_3.txt", package = "methylSig"),
              system.file("extdata", "AML_4.txt", package = "methylSig"), 
              system.file("extdata", "NBM_1.txt", package = "methylSig"),
              system.file("extdata", "NBM_2.txt", package = "methylSig"),
              system.file("extdata", "NBM_3.txt", package = "methylSig"),
              system.file("extdata", "NBM_4.txt", package = "methylSig"))


sample.id = c("AML1","AML2","AML3","AML4","NBM1","NBM2","NBM3","NBM4")

treatment = c(1,1,1,1,0,0,0,0)
#### Read Data ####
meth <- methylSigReadData(fileList, sample.ids = sample.id,
             assembly = "hg18", treatment = treatment, 
             context = "CpG", destranded=TRUE)
}
