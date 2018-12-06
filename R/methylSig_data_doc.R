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
#' @name CT_SNPs_hg19
#' @format A GenomicRanges object of length 1,321,463
#' @source \url{ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz}
NULL

#' Sample data
#'
#' The \code{sample_data} object contains the following items:
#' \describe{
#'   \item{meth}{ A \code{BSseq-class} object containing 6 samples total, with three in each group. Genome is hg19. }
#'   \item{tiled_meth}{ A tiled version of the \code{BSseq-class} object called \code{meth}. Tiles are 1000bp. Genome is hg19. }
#'   \item{msig_cpgs}{ A \code{GRanges-class} object containing the results of \code{methylSigCalc} on \code{data}. }
#'   \item{msig_tiles}{ A \code{GRanges-class} object containing the results of \code{methylSigCalc} on \code{tiled_meth}. }
#'   \item{tfbs}{ A \code{GRanges-class} object representing transcription factor binding sites. For use in \code{methylSig.tfbsEnrichTest}. Genome is hg19. }
#' }
#'
#' @docType data
#' @keywords datasets
#' @name sample_data
#' @format A mixture of BSseq and GRanges class objects. See documentation for test_data.
NULL
