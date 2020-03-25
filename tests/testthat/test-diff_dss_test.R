data(BS.cancer.ex, package = 'bsseqData')

bs = filter_loci_by_group_coverage(
    bs = BS.cancer.ex,
    group_column = 'Type',
    min_samples_per_group = c('cancer' = 2, 'normal' = 2))

pData(bs)$num_cov = c(9, 8, 10, 1, 3, 2)

small_test = bs[1:50]

bs_tile = tile_by_windows(bs, win_size = 5000)

bs_tile = filter_loci_by_group_coverage(
    bs = bs_tile,
    group_column = 'Type',
    min_samples_per_group = c('cancer' = 2, 'normal' = 2))

small_test_tile = bs_tile[1:50]

diff_fit = diff_dss_fit(
    bs = small_test,
    design = pData(small_test),
    formula = '~ Type')

#####################################

test_that('bs missing check', {
    expect_error(
        diff_dss_test(),
        'Must pass bs as a BSseq object',
        fixed = TRUE
    )
})

test_that('diff_fit missing check', {
    expect_error(
        diff_dss_test(bs = small_test),
        'Must pass diff_fit',
        fixed = TRUE
    )
})

test_that('contrast missing check', {
    expect_error(
        diff_dss_test(bs = small_test, diff_fit = diff_fit),
        'Must pass contrast',
        fixed = TRUE
    )
})

#####################################

test_that('bs type check', {
    expect_error(
        diff_dss_test(
            bs = 'blue',
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1)),
        'bs must be',
        fixed = TRUE
    )
})

test_that('diff_fit type check', {
    expect_error(
        diff_dss_test(
            bs = small_test,
            diff_fit = 'blue',
            contrast = matrix(c(0,1), ncol = 1)),
        'diff_fit must be a list.',
        fixed = TRUE
    )
})

test_that('diff_fit a list with correct names check', {
    expect_error(
        diff_dss_test(
            bs = small_test,
            diff_fit = list('a' = 'hello', 'b' = 'goodbye'),
            contrast = matrix(c(0,1), ncol = 1)),
        'diff_fit must be a list returned from diff_dss_fit',
        fixed = TRUE
    )
})

#####################################

test_that('Valid methylation_group_column name check', {
    expect_error(
        diff_dss_test(
            bs = small_test,
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_group_column = 'blue'),
        'not in column names of diff_fit$design',
        fixed = TRUE
    )
})

test_that('methylation_groups and methylation_group_column check', {
    expect_error(
        diff_dss_test(
            bs = small_test,
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_groups = c('case' = 'blue', 'control' = 'read')),
        'If methylation_groups is specified',
        fixed = TRUE
    )
})

test_that('methylation_groups type check', {
    expect_error(
        diff_dss_test(
            bs = small_test,
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_group_column = 'Type',
            methylation_groups = 2),
        'methylation_groups must be a named character vector',
        fixed = TRUE
    )
})

test_that('methylation_groups type check', {
    expect_error(
        diff_dss_test(
            bs = small_test,
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_group_column = 'Type',
            methylation_groups = c('blue' = 'blue', 'red' = 'red')),
        'methylation_groups must be a named vector with names',
        fixed = TRUE
    )
})

test_that('methylation_groups and methylation_group_column check', {
    expect_error(
        diff_dss_test(
            bs = small_test,
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_group_column = 'Type',
            methylation_groups = c('case' = 'blue', 'control' = 'red')),
        'Not all methylation_groups are in methylation_group_column',
        fixed = TRUE
    )
})

#####################################

test_that('Valid return, simple model, group methylation check', {
    diff_fit = diff_dss_fit(
        bs = small_test,
        design = pData(small_test),
        formula = '~ Type')

    diff_gr = diff_dss_test(
        bs = small_test,
        diff_fit = diff_fit,
        contrast = matrix(c(0,1), ncol = 1),
        methylation_group_column = 'Type',
        methylation_groups = c('case' = 'cancer', 'control' = 'normal')
    )

    expect_true(is(diff_gr, 'GRanges'))

})

test_that('Valid return, more complex model, no methylation check', {
    diff_fit = diff_dss_fit(
        bs = small_test,
        design = pData(small_test),
        formula = '~ Type + Pair')

    diff_gr = diff_dss_test(
        bs = small_test,
        diff_fit = diff_fit,
        contrast = matrix(c(0,1,0,0), ncol = 1)
    )

    expect_true(is(diff_gr, 'GRanges'))

})

test_that('Valid return, more complex model, methylation check', {
    diff_fit = diff_dss_fit(
        bs = small_test,
        design = pData(small_test),
        formula = '~ Type + num_cov')

    diff_gr = diff_dss_test(
        bs = small_test,
        diff_fit = diff_fit,
        contrast = matrix(c(0,1,0), ncol = 1),
        methylation_group_column = 'Type',
        methylation_groups = c('case' = 'cancer', 'control' = 'normal')
    )

    expect_true(is(diff_gr, 'GRanges'))

})

test_that('Valid return, numerical covariate model, percentile methylation check', {
    diff_fit = diff_dss_fit(
        bs = small_test,
        design = pData(small_test),
        formula = '~ num_cov')

    diff_gr = diff_dss_test(
        bs = small_test,
        diff_fit = diff_fit,
        contrast = matrix(c(0,1), ncol = 1),
        methylation_group_column = 'num_cov'
    )

    expect_true(is(diff_gr, 'GRanges'))

})

test_that('Valid return, simple model tiled, methylation check', {
    diff_fit = diff_dss_fit(
        bs = small_test_tile,
        design = pData(small_test_tile),
        formula = '~ Type')

    diff_gr = diff_dss_test(
        bs = small_test_tile,
        diff_fit = diff_fit,
        contrast = matrix(c(0,1), ncol = 1),
        methylation_group_column = 'Type',
        methylation_groups = c('case' = 'cancer', 'control' = 'normal')
    )

    expect_true(is(diff_gr, 'GRanges'))

})
