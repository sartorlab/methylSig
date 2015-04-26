# Testing of methylSig package
    library(methylSig)

    data(CT_SNPs_hg19)
    CT_SNPs_hg19_chr21 = CT_SNPs_hg19[which(seqnames(CT_SNPs_hg19) == 'chr21')]

# Relevant files are located at ~/latte/Methylation/Data/methylsig_testing
    setwd('~/latte/Methylation/Data/methylsig_testing')

# Files used to test methylSigReadData (only chr21)
#   IDH2_1_CpG.txt
#   IDH2_2_CpG.txt
#   IDH2_3_CpG.txt
#
#   NBM_1_CpG.txt
#   NBM_2_CpG.txt
#   NBM_3_CpG.txt
    msrd_files = list.files(pattern='CpG.txt.gz')

# Files used to test readBismarkData (only chr21)
#   IDH2mut_1_chr21_errbs.fastq_bismark_sorted.bismark.cov
#   IDH2mut_1_chr21_errbs.fastq_bismark_sorted.CpG_report.txt
#   IDH2mut_2_chr21_errbs.fastq_bismark_sorted.bismark.cov
#   IDH2mut_2_chr21_errbs.fastq_bismark_sorted.CpG_report.txt
#   IDH2mut_3_chr21_errbs.fastq_bismark_sorted.bismark.cov
#   IDH2mut_3_chr21_errbs.fastq_bismark_sorted.CpG_report.txt
#
#   NBM_1_chr21_errbs.fastq_bismark_sorted.bismark.cov
#   NBM_1_chr21_errbs.fastq_bismark_sorted.CpG_report.txt
#   NBM_2_chr21_errbs.fastq_bismark_sorted.bismark.cov
#   NBM_2_chr21_errbs.fastq_bismark_sorted.CpG_report.txt
#   NBM_3_chr21_errbs.fastq_bismark_sorted.bismark.cov
#   NBM_3_chr21_errbs.fastq_bismark_sorted.CpG_report.txt
    rbd_cov_files = list.files(pattern='bismark.cov.gz')
    rbd_cyt_files = list.files(pattern='CpG_report.txt')

#   Test for differences in choice of read function
    test_msrd = methylSigReadData(
        fileList = msrd_files,
        sample.ids = c('IDH2_1','IDH2_2','IDH2_3','NBM_1','NBM_2','NBM_3'),
        assembly = 'hg19',
        pipeline = 'bismark and methylKit',
        header = TRUE,
        context = 'CpG',
        resolution = "base",
        treatment = c(1,1,1,0,0,0),
        destranded = TRUE,
        maxCount = 500,
        minCount = 10,
        filterSNPs = FALSE,
        num.cores = 6,
        quiet = FALSE)
    test_rbd = readBismarkData(
        bismarkCovFiles = rbd_cov_files,
        cytosineCovFiles = rbd_cyt_files,
        sample.ids = c('IDH2_1','IDH2_2','IDH2_3','NBM_1','NBM_2','NBM_3'),
        assembly = 'hg19',
        pipeline = 'bismark',
        context = 'CpG',
        resolution = "base",
        treatment = c(1,1,1,0,0,0),
        destranded = TRUE,
        maxCount = 500,
        minCount = 10,
        filterSNPs = FALSE,
        num.cores = 6,
        quiet = FALSE)

    # Expectations
    # This is FALSE. Possibly because of phred quality step (ignore < 20) when
    # methylKit processes the SAM files. bismark_methylation_extractor does not
    # have an option like this. How to solve?
    nrow(test_msrd@data.coverage) == nrow(test_rbd@data.coverage)

# Testing binomialDiffCalc
#   Test different values of min.per.group
    test_msrd = methylSigReadData(
        fileList = msrd_files,
        sample.ids = c('IDH2_1','IDH2_2','IDH2_3','NBM_1','NBM_2','NBM_3'),
        assembly = 'hg19',
        pipeline = 'bismark and methylKit',
        header = TRUE,
        context = 'CpG',
        resolution = "base",
        treatment = c(1,1,1,0,0,0),
        destranded = TRUE,
        maxCount = 500,
        minCount = 10,
        filterSNPs = FALSE,
        num.cores = 6,
        quiet = FALSE)
    test_rbd = readBismarkData(
        bismarkCovFiles = rbd_cov_files,
        cytosineCovFiles = rbd_cyt_files,
        sample.ids = c('IDH2_1','IDH2_2','IDH2_3','NBM_1','NBM_2','NBM_3'),
        assembly = 'hg19',
        pipeline = 'bismark',
        context = 'CpG',
        resolution = "base",
        treatment = c(1,1,1,0,0,0),
        destranded = TRUE,
        maxCount = 500,
        minCount = 10,
        filterSNPs = FALSE,
        num.cores = 6,
        quiet = FALSE)

    test_bdf_msrd = binomialDiffCalc(
        meth = test_msrd,
        groups = c("Treatment"=1,"Control"=0),
        min.per.group = c(3,3))
    test_bdf_rbd = binomialDiffCalc(
        meth = test_rbd,
        groups = c("Treatment"=1,"Control"=0),
        min.per.group = c(3,3))

    # Expectations
    class(test_bdf_msrd)[1] == 'methylSigDiff'
    class(test_bdf_rbd)[1] == 'methylSigDiff'

# Testing methylSigCalc (resolution='site')
#   Test dispersion='Treatment', dispersion='Control', dispersion=1, dispersion=0, and dispersion='both'
#   Test local.disp=T and local.disp=F
#   When local.disp=T test different values of winsize.disp
#   Test local.meth=T and local.meth=F
#   When local.meth=T test different values of winsize.meth
#   Test different values of min.per.group
#   Test different weight functions
#   Test T.approx=T and T.approx=F
#   Test num.cores=1 and num.cores > 1
#   Test that an object of methylSigDiff-class is output
#
# Testing methylSigTile
#   Test different values of win.size
#   Test that an object of methylSigData-class is output with resolution=='region'
#
# Testing methylSigCalc (resolution='region')
#   Test dispersion='Treatment', dispersion='Control', dispersion=1, dispersion=0, and dispersion='both'
#   Test local.disp=T and local.disp=F
#   When local.disp=T test different values of winsize.disp
#   Test local.meth=T and local.meth=F
#   When local.meth=T test different values of winsize.meth
#   Test different values of min.per.group
#   Test different weight functions
#   Test T.approx=T and T.approx=F
#   Test num.cores=1 and num.cores > 1
#   Test that an object of methylSigDiff-class is output
#
# Testing getTFBSinfo
#   Not sure what to test
#
# Testing methylSig.tfbsEnrichTest
#
#
# Testing methylSigTileTFBS
#
