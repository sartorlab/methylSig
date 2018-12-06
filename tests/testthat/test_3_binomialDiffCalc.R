context('Test binomialDiffCalc')

utils::data(sample_data, package = 'methylSig')

test_that('Test CpG no local both dispersion', {
    result_binomial = binomialDiffCalc(
        meth = meth,
        comparison = 'DR_vs_DS',
        min.per.group = c(3,3))

    expect_true(is(result_binomial, 'GRanges'))
    expect_match(S4Vectors::metadata(result_binomial)$method, 'binomialDiffCalc')
})
