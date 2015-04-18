\name{methylSig-package}
\alias{methylSig-package}
\alias{methylSig}
\docType{package}
\title{
a whole genome DNA methylation analysis pipline
}
\description{
MethylSig is a method for testing differential methylated cytosines (DMCs) or regions (DMRs) in whole-genome bisulfite sequencing (bis-seq) or reduced representation bisulfite sequencing (RRBS) experiments. MethylSig uses a beta binomial model to test for significant differences between groups of samples. Several options exist for either site-specific or sliding window tests, combining strands, and for variance estimation. It allows annotating the resulting regions to multiple genome features, and visualizing the results for chosen genomic regions.
}
\details{
\tabular{ll}{
Package: \tab methylSig\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2014-04-05\cr
License: \tab GPL-3\cr
Depends: \tab methods, boot, parallel \cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Yongseok Park, Maureen A. Sartor.

Maintainer: Who to complain to <yongpark@pitt.edu>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
\examples{

library(methylSig)

##############  MethylSig ################
############### Example ##################

cpgInfo = getCpGInfo(system.file("annotation", "cpgi.hg18.bed.txt", package = "methylSig"))
tfbsInfo = getTFBSInfo(system.file("annotation", "tfbsUniform.txt", package = "methylSig"))
refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt", package = "methylSig"))

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
meth <- methylSigReadData(fileList, sample.ids = sample.id,assembly = "hg18", 
             treatment = treatment, context = "CpG", destranded=TRUE,  num.cores = 5)

#### Differential methylation analysis #####
myDiffSigboth = methylSigCalc(meth, groups=c(1,0), min.per.group=3, num.cores=5)

####
###### Variance from sample treatment group "0" only
myDiffSignorm = methylSigCalc(meth, groups=c(1,0), dispersion=0, min.per.group=3, num.cores=5)

####
###### Variance from sample treatment group "0" only and local information
myDiffSignormLoc = methylSigCalc(meth, groups=c(1,0), dispersion=0, min.per.group=3, local.disp=TRUE, num.cores=5)

### Tiled analysis
methTile = methylSigTile(meth)

myDiffSigbothTile = methylSigCalc(methTile, groups=c(1,0), min.per.group=3, num.cores=5)
myDiffSignormTileLoc = methylSigCalc(methTile, groups=c(1,0), dispersion=0, min.per.group=3, local.disp=TRUE, num.cores=5)

#################
##   CpG plots
#################
cpgAnn = cpgAnnotation(cpgInfo,myDiffSigboth)
cpgAnnDmc = cpgAnnotation(cpgInfo,myDiffSigboth[myDiffSigboth[,"qvalue"] < 0.05 
                                            & abs(myDiffSigboth[,"meth.diff"])>25])
cpgAnnotationPlot(cpgAnn,main="ALL")
cpgAnnotationPlot(cpgAnnDmc, main="DMC")

refGeneAnn = refGeneAnnotation(refGeneInfo,myDiffSigboth)
refGeneAnnDmc = refGeneAnnotation(refGeneInfo,myDiffSigboth[myDiffSigboth[,"qvalue"] < 0.05 
                                            & abs(myDiffSigboth[,"meth.diff"])>25])


refGeneAnnotationPlot(refGeneAnn,main="ALL",priority=c("promoter","cds", "noncoding", "5'utr", "3'utr"))
refGeneAnnotationPlot(refGeneAnnDmc, main="DMC",priority=c("promoter","cds", "noncoding", "5'utr", "3'utr"))

############ tfbs Enrichment Test ###############
pvalue = methylSig.tfbsEnrichTest(myDiffSigboth,
                   myDiffSigboth[,"qvalue"] < 0.05 & abs(myDiffSigboth[,"meth.diff"]) > 25, tfbsInfo)  

############ plot data ####################
methylSigPlot(meth, "chr21", c(43000000, 44000000), groups=c(1,0), cpgInfo=cpgInfo,refGeneInfo=refGeneInfo,
    myDiff=myDiffSigboth,tfbsInfo=tfbsInfo,tfbsDense=FALSE,sigQ=0.05)


methylSigPlot(meth, "chr21", c(43800000, 43900000), groups=c(1,0), cpgInfo=cpgInfo,refGeneInfo=refGeneInfo,
    myDiff=myDiffSigboth,tfbsInfo=tfbsInfo,tfbsDense=FALSE,sigQ=0.05)

}

\keyword{internal}
