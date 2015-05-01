context('Test methylSigReadData')

# CT_SNPs for chr21 only
    data(CT_SNPs_hg19)
    CT_SNPs_hg19_chr21 = CT_SNPs_hg19[seqnames(CT_SNPs_hg19) == 'chr21']

# Files used to test methylSigReadData (only chr21)
#   IDH2_1_CpG.txt.gz
#   IDH2_2_CpG.txt.gz
#   IDH2_3_CpG.txt.gz
#
#   NBM_1_CpG.txt.gz
#   NBM_2_CpG.txt.gz
#   NBM_3_CpG.txt.gz
    msrd_files = c(
        system.file('extdata', 'IDH2_1_CpG.txt.gz', package='methylSig'),
        system.file('extdata', 'IDH2_2_CpG.txt.gz', package='methylSig'),
        system.file('extdata', 'IDH2_3_CpG.txt.gz', package='methylSig'),
        system.file('extdata', 'NBM_1_CpG.txt.gz', package='methylSig'),
        system.file('extdata', 'NBM_2_CpG.txt.gz', package='methylSig'),
        system.file('extdata', 'NBM_3_CpG.txt.gz', package='methylSig')
        )

# Testing methylSigReadData
#   Test different values of minCount and maxCount
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
        minCount = 20,
        filterSNPs = FALSE,
        num.cores = 2,
        quiet = FALSE)

    # Expectations
    # The minimum coverage test is actually dependent on previous steps
    # for instance bismark_methylation_extractor and methylKit both have
    # minimum coverage parameters. For this test data, the cutoff was 10
    # in both cases, meaning the minimum will always be 10 when the
    # minCount <= 10.
    test_that('Filtering by minCount works',{
        expect_more_than(min(test_msrd@data.coverage,na.rm=T), 19)
    })

    # The maximum coverage test is for 2*maxCount because destranding happens
    # AFTER the maxCount filter. Meaning sites can have up to 2*maxCount reads.
    test_that('Filtering by maxCoutn works',{
        expect_less_than(max(test_msrd@data.coverage,na.rm=T), 500*2 + 1)
    })

    test_that('Object of methylSigData class is returned',{
        expect_equal(class(test_msrd)[1], 'methylSigData')
    })

#   Test destranded=F
    test_msrd = methylSigReadData(
        fileList = msrd_files,
        sample.ids = c('IDH2_1','IDH2_2','IDH2_3','NBM_1','NBM_2','NBM_3'),
        assembly = 'hg19',
        pipeline = 'bismark and methylKit',
        header = TRUE,
        context = 'CpG',
        resolution = "base",
        treatment = c(1,1,1,0,0,0),
        destranded = FALSE,
        maxCount = 500,
        minCount = 10,
        filterSNPs = FALSE,
        num.cores = 2,
        quiet = FALSE)

    # Expectations
    test_that('Filtering by minCount when not destranded works',{
        expect_equal(min(test_msrd@data.coverage,na.rm=T), 10)
    })

    # In the case of destranded=F, the max should really be max (instead of
    # 2*max) because destranding occurs after maxCount filtering. As noted.
    test_that('Filtering by maxCount when not destranded works',{
        expect_less_than(max(test_msrd@data.coverage,na.rm=T), 501)
    })

    test_that('Object of methylSigData class is returned when destranded',{
        expect_equal(class(test_msrd)[1], 'methylSigData')
    })

#   Test filterSNPs=T and destranded=F
    test_msrd = methylSigReadData(
        fileList = msrd_files,
        sample.ids = c('IDH2_1','IDH2_2','IDH2_3','NBM_1','NBM_2','NBM_3'),
        assembly = 'hg19',
        pipeline = 'bismark and methylKit',
        header = TRUE,
        context = 'CpG',
        resolution = "base",
        treatment = c(1,1,1,0,0,0),
        destranded = FALSE,
        maxCount = 500,
        minCount = 10,
        filterSNPs = TRUE,
        num.cores = 2,
        quiet = FALSE)

    # Expectations
    test_that('filterSNPs works when not destranded',{
        expect_equal(any(test_msrd@data.start %in% start(CT_SNPs_hg19_chr21)), FALSE)
    })

    test_that('Object of methylSigData class is returned when not destranded and filterSNPs',{
        expect_equal(class(test_msrd)[1], 'methylSigData')
    })

#   Test filterSNPs=T and destranded=T
#       The any(test_msrd@data.start %in% start(CT_SNPs_hg19_chr21)) == FALSE
#       test fails because destranding shifts C > T + 1 sites back to the C > T
#       site.
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
        filterSNPs = TRUE,
        num.cores = 2,
        quiet = FALSE)

    # Expectations
    test_that('filterSNPs works when destranded',{
        expect_equal(any(test_msrd@data.start %in% start(CT_SNPs_hg19_chr21)), FALSE)
    })

    test_that('Object of methylSigData class is returned destranded and filterSNPs',{
        expect_equal(class(test_msrd)[1], 'methylSigData')
    })
