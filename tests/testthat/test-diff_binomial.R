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
