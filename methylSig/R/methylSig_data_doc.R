#' CpG Index for hg19
#'
#' A GenomicsRanges object giving the coordinates (in hg19) of all C > T SNPs.
#' Start coordinates are 0-based and end coordinates are 1-based. Starting from
#' 1000 Genomes Data we used \code{bcftools filter} with \code{-i 'AF[0]>0.05'}
#' to pull all sites with alternate frequency greater than 0.05. We then used
#' \code{grep -P '(C\tT)'} and \code{grep -P '(VT=SNP)'} to collect all C > T
#' SNPs.
#'
#' @docType data
#' @keywords datasets
#' @name data_CT_SNPs_hg19
#' @aliases CT_SNPs_hg19
#' @usage data(CT_SNPs_hg19)
#' @format A GenomicRanges object of length 1,321,463
#' @source \url{ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz}
NULL

#' Example methylSigData and methylSigDiff objects
#'
#' \code{meth} is a \code{\link{methylSigData}} object containing 8 samples total.
#' There are 4 AML (acute myeloid leukemia) samples and 4 NBM (normal bone marrow)
#' samples. Data are in hg18 coordinates.
#' \code{mySigDiffboth} is a \code{\link{methylSigDiff}} object storing the result
#' of the site-wise differential methylation test using \code{\link{methylSigCalc}}.
#'
#' @docType data
#' @keywords datasets
#' @name sampleData
#' @usage data(sampleData)
NULL

#' Example methylSigData object
#'
#' \code{meth} is a \code{\link{methylSigData}} object containing 8 samples total.
#' There are 4 AML (acute myeloid leukemia) samples and 4 NBM (normal bone marrow)
#' samples. Data are in hg18 coordinates.
#'
#' @docType data
#' @keywords datasets
#' @name meth
#' @usage data(sampleData)
NULL

#' Example methylSigDiff object
#'
#' \code{mySigDiffboth} is a \code{\link{methylSigDiff}} object storing the result
#' of the site-wise differential methylation test using \code{\link{methylSigCalc}}.
#'
#' @docType data
#' @keywords datasets
#' @name myDiffSigboth
#' @usage data(sampleData)
NULL
