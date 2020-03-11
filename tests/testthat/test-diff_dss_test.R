data(BS.cancer.ex, package = 'bsseqData')

bs = filter_loci_by_group_coverage(
    bs = BS.cancer.ex,
    group_column = 'Type',
    c('cancer' = 2, 'normal' = 2))

small_test = bs[1:50]

diff_fit = diff_dss_fit(
    bs = small_test,
    design = pData(bs),
    formula = '~ Type')

#####################################

test_that('diff_fit missing check', {
    expect_error(
        diff_dss_test(),
        'Must pass diff_fit',
        fixed = TRUE
    )
})

test_that('contrast missing check', {
    expect_error(
        diff_dss_test(diff_fit = diff_fit),
        'Must pass contrast',
        fixed = TRUE
    )
})

#####################################

test_that('diff_fit a list check', {
    expect_error(
        diff_dss_test(
            diff_fit = 'blue',
            contrast = matrix(c(0,1), ncol = 1)),
        'diff_fit must be a list.',
        fixed = TRUE
    )
})

test_that('diff_fit a list with correct names check', {
    expect_error(
        diff_dss_test(
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
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_group_column = 'blue'),
        'not in column names of diff_fit$design', # () seem to be a problem
        fixed = TRUE
    )
})

test_that('methylation_groups and methylation_group_column check', {
    expect_error(
        diff_dss_test(
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_groups = c('case' = 'blue', 'control' = 'read')),
        'If methylation_groups is specified', # () seem to be a problem
        fixed = TRUE
    )
})

test_that('methylation_groups type check', {
    expect_error(
        diff_dss_test(
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_group_column = 'Type',
            methylation_groups = 2),
        'methylation_groups must be a named character vector', # () seem to be a problem
        fixed = TRUE
    )
})

test_that('methylation_groups type check', {
    expect_error(
        diff_dss_test(
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_group_column = 'Type',
            methylation_groups = c('blue' = 'blue', 'red' = 'red')),
        'methylation_groups must be a named vector with names', # () seem to be a problem
        fixed = TRUE
    )
})

test_that('methylation_groups and methylation_group_column check', {
    expect_error(
        diff_dss_test(
            diff_fit = diff_fit,
            contrast = matrix(c(0,1), ncol = 1),
            methylation_group_column = 'Type',
            methylation_groups = c('case' = 'blue', 'control' = 'red')),
        'Not all methylation_groups are in methylation_group_column', # () seem to be a problem
        fixed = TRUE
    )
})

#####################################

# result = diff_dss_test(
#     diff_fit = diff_fit,
#     contrast = matrix(c(0,1), ncol = 1)
# )
