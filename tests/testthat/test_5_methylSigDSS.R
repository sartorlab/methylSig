context('Test methylSigDSS')

utils::data(sample_data, package = 'methylSig')

design1 = data.frame(
    group = bsseq::pData(meth)$DR_vs_DS)

design2 = data.frame(
    group = bsseq::pData(meth)$DR_vs_DS,
    subject = factor(c(1,1,2,2,3,3)))

test_that('Test with intercept', {
    contrast = matrix(c(0,1), ncol = 1)
    result_dss = methylSigDSS(
        meth = meth,
        design = design1,
        formula = '~ group',
        contrast = contrast,
        group.term = 'group',
        min.per.group=c(3,3))

    expect_match(class(result_dss), 'GRanges')
    expect_match(S4Vectors::metadata(result_dss)$method, 'methylSigDSS')
})

test_that('Test without intercept', {
    contrast = matrix(c(-1,1), ncol = 1)
    result_dss = methylSigDSS(
        meth = meth,
        design = design1,
        formula = '~ 0 + group',
        contrast = contrast,
        group.term = 'group',
        min.per.group=c(3,3))

    expect_match(class(result_dss), 'GRanges')
    expect_match(S4Vectors::metadata(result_dss)$method, 'methylSigDSS')
})

test_that('Test similar to first but with extra design columns', {
    contrast = matrix(c(0,1), ncol = 1)
    result_dss = methylSigDSS(
        meth = meth,
        design = design2,
        formula = '~ group',
        contrast = contrast,
        group.term = 'group',
        min.per.group=c(3,3))

    expect_match(class(result_dss), 'GRanges')
    expect_match(S4Vectors::metadata(result_dss)$method, 'methylSigDSS')
})

test_that('Test multiple formula terms', {
    contrast = matrix(c(0,1,0,0), ncol = 1)
    result_dss = methylSigDSS(
        meth = meth,
        design = design2,
        formula = '~ group + subject',
        contrast = contrast,
        group.term = 'group',
        min.per.group=c(3,3))

    expect_match(class(result_dss), 'GRanges')
    expect_true('meth.2' %in% colnames(GenomicRanges::mcols(result_dss)))
    expect_match(S4Vectors::metadata(result_dss)$method, 'methylSigDSS')
})

test_that('Test alternate formula', {
    contrast = matrix(c(0,0,0,1), ncol = 1)
    result_dss = methylSigDSS(
        meth = meth,
        design = design2,
        formula = '~ subject + group',
        contrast = contrast,
        group.term = 'group',
        min.per.group=c(3,3))

    expect_match(class(result_dss), 'GRanges')
    expect_match(S4Vectors::metadata(result_dss)$method, 'methylSigDSS')
})
