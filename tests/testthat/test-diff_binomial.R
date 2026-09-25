data(BS.cancer.ex, package = 'bsseqData')

bs = filter_loci_by_group_coverage(
    bs = BS.cancer.ex,
    group_column = 'Type',
    c('cancer' = 2, 'normal' = 2))

small_test = bs[1:50]

small_test_tile = tile_by_windows(bs = small_test, win_size = 5000)

#####################################

test_that('bs missing check', {
    expect_error(
        diff_binomial(),
        'Must pass bs as a BSseq object',
        fixed = TRUE
    )
})

test_that('group_column missing check', {
    expect_error(
        diff_binomial(bs = small_test),
        'Must pass group_column',
        fixed = TRUE
    )
})

test_that('comparison_groups missing check', {
    expect_error(
        diff_binomial(
            bs = small_test,
            group_column = 'Type'),
        'Must pass comparison_groups',
        fixed = TRUE
    )
})

#####################################

test_that('bs type check', {
    expect_error(
        diff_binomial(
            bs = 'blue',
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal')),
        'bs must be',
        fixed = TRUE
    )
})

test_that('group_column type check', {
    expect_error(
        diff_binomial(
            bs = small_test,
            group_column = c(1, 3),
            comparison_groups = c('case' = 'cancer', 'control' = 'normal')),
        'group_column must be',
        fixed = TRUE
    )
})

test_that('comparison_groups type check', {
    expect_error(
        diff_binomial(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 1, 'control' = 2)),
        'comparison_groups must be',
        fixed = TRUE
    )
})

#####################################

test_that('Valid group_column name check', {
    expect_error(
        diff_binomial(
            bs = small_test,
            group_column = 'blue',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal')),
        'not in column names of pData(bs):',
        fixed = TRUE
    )
})

test_that('Valid comparison_groups values check', {
    expect_error(
        diff_binomial(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'blue', 'control' = 'normal')),
        'Not all comparison_groups are in group_column',
        fixed = TRUE
    )
})

test_that('Valid comparison_groups name check', {
    expect_error(
        diff_binomial(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('blue' = 'cancer', 'control' = 'normal')),
        'comparison_groups vector must be a named vector with',
        fixed = TRUE
    )
})

#####################################

test_that('Test 1', {
    diff_gr = diff_binomial(
        bs = small_test,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'))

    expect_true(is(diff_gr, 'GRanges'))
})

test_that('Test 2', {
    diff_gr = diff_binomial(
        bs = small_test_tile,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'))

    expect_true(is(diff_gr, 'GRanges'))
})

test_that('One sample in each group', {
    one_each = small_test[, c(1, 4)]
    cov = as.matrix(getCoverage(one_each, type = 'Cov'))
    n_untestable = sum(cov[, 1] == 0 | cov[, 2] == 0)

    diff_gr = suppressMessages(diff_binomial(
        bs = one_each,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal')))

    expect_true(is(diff_gr, 'GRanges'))
    expect_equal(length(diff_gr), length(small_test) - n_untestable)
})

test_that('Samples outside case and control are not used', {
    three_groups = small_test
    pData(three_groups)$Group = c('cancer', 'cancer', 'other', 'normal', 'normal', 'other')
    comparison_groups = c('case' = 'cancer', 'control' = 'normal')

    with_other = suppressMessages(diff_binomial(three_groups, 'Group', comparison_groups))
    without_other = suppressMessages(diff_binomial(
        three_groups[, pData(three_groups)$Group != 'other'], 'Group', comparison_groups))

    expect_equal(with_other, without_other)
})

test_that('log_lik_ratio is the binomial deviance difference', {
    diff_gr = suppressMessages(diff_binomial(
        bs = small_test,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal')))

    cov = as.matrix(getCoverage(small_test, type = 'Cov'))
    meth = as.matrix(getCoverage(small_test, type = 'M'))
    group = factor(pData(small_test)$Type)
    idx = match(start(diff_gr), start(small_test))

    deviance_diff = vapply(idx, function(i) {
        reads = cbind(meth[i, ], cov[i, ] - meth[i, ])
        keep = cov[i, ] > 0
        # glm() warns when a group's fitted proportion is 0 or 1
        suppressWarnings(
            stats::deviance(stats::glm(reads[keep, ] ~ 1, family = stats::binomial)) -
            stats::deviance(stats::glm(reads[keep, ] ~ group[keep], family = stats::binomial)))
    }, 1)

    expect_equal(diff_gr$log_lik_ratio, deviance_diff, tolerance = 1e-6)
})

test_that('meth_diff is computed before rounding', {
    diff_gr = suppressMessages(diff_binomial(
        bs = small_test,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal')))

    cov = as.matrix(getCoverage(small_test, type = 'Cov'))
    meth = as.matrix(getCoverage(small_test, type = 'M'))
    case = pData(small_test)$Type == 'cancer'
    idx = match(start(diff_gr), start(small_test))
    meth_case = rowSums(meth[idx, case]) / rowSums(cov[idx, case]) * 100
    meth_control = rowSums(meth[idx, !case]) / rowSums(cov[idx, !case]) * 100

    expect_equal(diff_gr$meth_diff, round(meth_case - meth_control, 2))
})

test_that('Loci without coverage in case or control are dropped', {
    no_case = small_test
    cov = as.matrix(getCoverage(no_case, type = 'Cov'))
    meth = as.matrix(getCoverage(no_case, type = 'M'))
    case = pData(no_case)$Type == 'cancer'
    cov[1:3, case] = 0
    meth[1:3, case] = 0
    no_case = BSseq(Cov = cov, M = meth, gr = granges(no_case), pData = pData(no_case),
        sampleNames = sampleNames(no_case))

    expect_message(
        diff_gr <- diff_binomial(
            bs = no_case,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal')),
        'loci were dropped because case or control has no coverage.',
        fixed = TRUE)
    expect_false(any(start(no_case)[1:3] %in% start(diff_gr)))
    expect_false(anyNA(diff_gr$pvalue))
    expect_false(anyNA(diff_gr$fdr))
})
