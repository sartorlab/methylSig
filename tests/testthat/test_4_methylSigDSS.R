context('Test methylSigDSS')

data(data, package = 'methylSig')

design1 = data.frame(
    group = bsseq::pData(data)$DR_vs_DS)

design2 = data.frame(
    group = bsseq::pData(data)$DR_vs_DS,
    subject = factor(c(1,1,2,2,3,3)))

# Test with intercept
contrast = matrix(c(0,1), ncol = 1)
result1 = methylSigDSS(
    meth = data,
    design = design1,
    formula = '~ group',
    contrast = contrast,
    group.term = 'group',
    min.per.group=c(3,3))

test_that('Returns GRanges', {
    expect_match(class(result1), 'GRanges')
})

# Test without intercept
contrast = matrix(c(-1,1), ncol = 1)
result2 = methylSigDSS(
    meth = data,
    design = design1,
    formula = '~ 0 + group',
    contrast = contrast,
    group.term = 'group',
    min.per.group=c(3,3))

test_that('Returns GRanges', {
    expect_match(class(result2), 'GRanges')
})

# Test that they are the same
test_that('Intercept and no-intercept tests are equivalent', {
    expect_equal(result1$pvalue, result2$pvalue, tolerance = 0.0000001)
})

contrast = matrix(c(0,1), ncol = 1)
result3 = methylSigDSS(
    meth = data,
    design = design2,
    formula = '~ group',
    contrast = contrast,
    group.term = 'group',
    min.per.group=c(3,3))

test_that('Returns GRanges', {
    expect_match(class(result3), 'GRanges')
})

contrast = matrix(c(0,1,0,0), ncol = 1)
result4 = methylSigDSS(
    meth = data,
    design = design2,
    formula = '~ group + subject',
    contrast = contrast,
    group.term = 'group',
    min.per.group=c(3,3))

test_that('Returns GRanges', {
    expect_match(class(result4), 'GRanges')
})

contrast = matrix(c(-1,1,0,0), ncol = 1)
result5 = methylSigDSS(
    meth = data,
    design = design2,
    formula = '~ 0 + group + subject',
    contrast = contrast,
    group.term = 'group',
    min.per.group=c(3,3))

test_that('Returns GRanges', {
    expect_match(class(result5), 'GRanges')
})

contrast = matrix(c(0,0,0,1), ncol = 1)
result6 = methylSigDSS(
    meth = data,
    design = design2,
    formula = '~ subject + group',
    contrast = contrast,
    group.term = 'group',
    min.per.group=c(3,3))

test_that('Returns GRanges', {
    expect_match(class(result6), 'GRanges')
})
