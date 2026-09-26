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

test_that('covariate_percentiles check', {
    diff_fit = diff_dss_fit(
        bs = small_test,
        design = pData(small_test),
        formula = '~ num_cov')

    for (bad in list(c(75, 25), 25, c(-1, 75), c(25, 101), c(25, NA), c('25', '75'))) {
        expect_error(
            diff_dss_test(
                bs = small_test,
                diff_fit = diff_fit,
                contrast = matrix(c(0,1), ncol = 1),
                methylation_group_column = 'num_cov',
                covariate_percentiles = bad),
            'covariate_percentiles must be two increasing numbers from 0 to 100.',
            fixed = TRUE
        )
    }
})

test_that('covariate_percentiles groups samples', {
    diff_fit = diff_dss_fit(
        bs = small_test,
        design = pData(small_test),
        formula = '~ num_cov')

    meth_mat = as.matrix(bsseq::getCoverage(small_test, type = 'M'))
    cov_mat = as.matrix(bsseq::getCoverage(small_test, type = 'Cov'))
    group_meth = function(idx) {
        round(rowSums(meth_mat[, idx, drop = FALSE]) / rowSums(cov_mat[, idx, drop = FALSE]) * 100, 2)
    }

    # num_cov is c(9, 8, 10, 1, 3, 2)
    # Default 25 and 75: case is num_cov <= 2.25 (1, 2), control is >= 8.75 (9, 10)
    default_gr = diff_dss_test(
        bs = small_test,
        diff_fit = diff_fit,
        contrast = matrix(c(0,1), ncol = 1),
        methylation_group_column = 'num_cov')
    expect_equal(unname(default_gr$meth_case), unname(group_meth(c(4, 6))))
    expect_equal(unname(default_gr$meth_control), unname(group_meth(c(1, 3))))

    # 50 and 50: case is num_cov <= 5.5 (1, 2, 3), control is >= 5.5 (8, 9, 10)
    median_gr = diff_dss_test(
        bs = small_test,
        diff_fit = diff_fit,
        contrast = matrix(c(0,1), ncol = 1),
        methylation_group_column = 'num_cov',
        covariate_percentiles = c(49.9, 50))
    expect_equal(unname(median_gr$meth_case), unname(group_meth(c(4, 5, 6))))
    expect_equal(unname(median_gr$meth_control), unname(group_meth(c(1, 2, 3))))

    # The test statistics don't depend on the grouping
    expect_equal(median_gr$stat, default_gr$stat)
})

test_that('covariate_percentiles that put samples in both groups', {
    diff_fit = diff_dss_fit(
        bs = small_test,
        design = pData(small_test),
        formula = '~ num_cov')
    # The 10 and 50 percentiles are both 1
    diff_fit$design$num_cov = c(1, 1, 1, 1, 1, 2)

    expect_error(
        diff_dss_test(
            bs = small_test,
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_group_column = 'num_cov',
            covariate_percentiles = c(10, 50)),
        'covariate_percentiles 10 and 50 of methylation_group_column num_cov put some samples in both groups.',
        fixed = TRUE
    )
})

test_that('contrast of the wrong size lists the columns of diff_fit$X', {
    for (contrast in list(c(0, 1, 0), matrix(1, ncol = 1))) {
        expect_error(
            diff_dss_test(bs = small_test, diff_fit = diff_fit, contrast = contrast),
            'contrast needs 2 rows (or a vector of length 2), one per column of diff_fit$X, in this order:\n  1: (Intercept)\n  2: Typenormal',
            fixed = TRUE
        )
    }
})

test_that('contrast type and value checks', {
    expect_error(
        diff_dss_test(bs = small_test, diff_fit = diff_fit, contrast = c('0', '1')),
        'contrast must be a numeric vector or matrix.',
        fixed = TRUE
    )
    expect_error(
        diff_dss_test(bs = small_test, diff_fit = diff_fit, contrast = c(0, NA)),
        'contrast must not have NA values.',
        fixed = TRUE
    )
    expect_error(
        diff_dss_test(bs = small_test, diff_fit = diff_fit, contrast = c(0, 0)),
        'Each column of contrast must have a nonzero value.',
        fixed = TRUE
    )
    expect_error(
        diff_dss_test(bs = small_test, diff_fit = diff_fit, contrast = cbind(c(0, 1), c(0, 2))),
        'The columns of contrast must be linearly independent.',
        fixed = TRUE
    )
    expect_error(
        diff_dss_test(bs = small_test, diff_fit = diff_fit, contrast = c('(Intercept)' = 0, 'Typecancer' = 1)),
        'The names of contrast ((Intercept), Typecancer) are not the columns of diff_fit$X.',
        fixed = TRUE
    )
})

test_that('contrast as a vector, a matrix, or by name gives the same test', {
    run = function(contrast) {
        suppressMessages(diff_dss_test(bs = small_test, diff_fit = diff_fit, contrast = contrast))
    }
    by_matrix = run(matrix(c(0, 1), ncol = 1))

    expect_equal(run(c(0, 1)), by_matrix)
    expect_equal(run(c('Typenormal' = 1, '(Intercept)' = 0)), by_matrix)
    expect_equal(run(matrix(c(1, 0), ncol = 1, dimnames = list(c('Typenormal', '(Intercept)'), NULL))), by_matrix)
})

test_that('contrast message says what is tested', {
    expect_message(
        diff_dss_test(bs = small_test, diff_fit = diff_fit, contrast = c(0, 1)),
        'Testing Typenormal = 0',
        fixed = TRUE
    )
    expect_message(
        diff_dss_test(bs = small_test, diff_fit = diff_fit, contrast = cbind(c(0, 1), c(1, -0.5))),
        'Testing Typenormal = 0 and (Intercept) - 0.5 * Typenormal = 0',
        fixed = TRUE
    )
})

test_that('A character methylation_group_column needs methylation_groups', {
    expect_error(
        diff_dss_test(
            bs = small_test,
            diff_fit = diff_fit,
            contrast = c(0, 1),
            methylation_group_column = 'Type'),
        'methylation_group_column Type is a character column, so methylation_groups must give its case and control values',
        fixed = TRUE
    )
})

test_that('methylation_group_column must be character, factor, or numeric', {
    logical_fit = diff_fit
    logical_fit$design$is_cancer = logical_fit$design$Type == 'cancer'

    expect_error(
        diff_dss_test(
            bs = small_test,
            diff_fit = logical_fit,
            contrast = c(0, 1),
            methylation_group_column = 'is_cancer'),
        'methylation_group_column is_cancer must be a character, factor, or numeric column of diff_fit$design, not logical.',
        fixed = TRUE
    )
})

test_that('methylation_groups is ignored for a numeric methylation_group_column', {
    num_fit = suppressMessages(diff_dss_fit(
        bs = small_test,
        design = pData(small_test),
        formula = '~ num_cov'))
    run = function(...) {
        suppressMessages(diff_dss_test(
            bs = small_test,
            diff_fit = num_fit,
            contrast = c(0, 1),
            methylation_group_column = 'num_cov',
            ...))
    }

    expect_warning(
        with_groups <- run(methylation_groups = c('case' = 'cancer', 'control' = 'normal')),
        'methylation_groups is ignored because methylation_group_column num_cov is numeric.',
        fixed = TRUE
    )
    expect_equal(with_groups, run())
})

test_that('Methylation rates come from the fit loci of bs, in order', {
    run = function(bs) {
        suppressMessages(diff_dss_test(
            bs = bs,
            diff_fit = diff_fit,
            contrast = c(0, 1),
            methylation_group_column = 'Type',
            methylation_groups = c('case' = 'cancer', 'control' = 'normal')))
    }
    expected = run(small_test)

    # More loci than were fit, and the fit loci in another order
    expect_equal(run(bs[1:100]), expected)
    expect_equal(run(small_test[rev(seq_along(small_test))]), expected)

    expect_error(
        run(small_test[1:40]),
        '10 of the 50 loci in diff_fit$gr are not in bs. Use the bs given to diff_dss_fit().',
        fixed = TRUE
    )
    expect_error(
        run(small_test[, 1:4]),
        'bs has 4 samples, but diff_fit$design has 6 rows. Use the bs given to diff_dss_fit().',
        fixed = TRUE
    )
})
