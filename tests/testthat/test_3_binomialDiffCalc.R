context('Test binomialDiffCalc')

data(data, package = 'methylSig')

test_that('Test CpG no local both dispersion', {
    result_binomial = binomialDiffCalc(
        meth = data,
        comparison = 'DR_vs_DS',
        min.per.group = c(3,3))

    expect_true(is(result_binomial, 'GRanges'))
    expect_match(S4Vectors::metadata(result_binomial)$method, 'binomialDiffCalc')
})
