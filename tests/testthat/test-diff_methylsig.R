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
        diff_methylsig(),
        'Must pass bs as a BSseq object',
        fixed = TRUE
    )
})

test_that('group_column missing check', {
    expect_error(
        diff_methylsig(bs = small_test),
        'Must pass group_column',
        fixed = TRUE
    )
})

test_that('comparison_groups missing check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type'),
        'Must pass comparison_groups',
        fixed = TRUE
    )
})

test_that('disp_groups missing check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal')),
        'Must pass disp_groups',
        fixed = TRUE
    )
})

#####################################

test_that('bs type check', {
    expect_error(
        diff_methylsig(
            bs = 'blue',
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'bs must be',
        fixed = TRUE
    )
})

test_that('group_column type check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = c(1, 3),
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'group_column must be',
        fixed = TRUE
    )
})

test_that('comparison_groups type check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 1, 'control' = 2),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'comparison_groups must be',
        fixed = TRUE
    )
})

test_that('disp_groups type check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = 1, 'control' = 2),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'disp_groups must be',
        fixed = TRUE
    )
})

test_that('local_window_size type check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 'a',
            t_approx = TRUE,
            n_cores = 1),
        'local_window_size must be',
        fixed = TRUE
    )
})

test_that('local_weight_function type check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            local_weight_function = 'b',
            t_approx = TRUE,
            n_cores = 1),
        'local_weight_function must be',
        fixed = TRUE
    )
})

test_that('t_approx type check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            t_approx = 'c',
            n_cores = 1),
        't_approx must be',
        fixed = TRUE
    )
})

test_that('n_cores type check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 'a'),
        'n_cores must be',
        fixed = TRUE
    )
})

#####################################

test_that('Valid group_column name check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'blue',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'not in column names of pData(bs):',
        fixed = TRUE
    )
})

test_that('Valid comparison_groups values check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'blue', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'Not all comparison_groups are in group_column',
        fixed = TRUE
    )
})

test_that('Valid comparison_groups name check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('blue' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'comparison_groups vector must be a named vector with',
        fixed = TRUE
    )
})

test_that('Valid disp_groups values check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = FALSE, 'control' = FALSE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'disp_groups must be a named logical vector with at least one TRUE value corresponding',
        fixed = TRUE
    )
})

test_that('Valid disp_groups name check', {
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('blue' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'disp_groups vector must be a named vector with names',
        fixed = TRUE
    )
})

test_that('Check for invalid local_window_size == 0 && regions state', {
    expect_error(
        diff_methylsig(
            bs = small_test_tile,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 50,
            t_approx = TRUE,
            n_cores = 1),
        'Cannot use local information on region-resolution data. Detected local_window_size',
        fixed = TRUE
    )
})

#####################################

test_that('Test 1', {
    diff_gr = diff_methylsig(
        bs = small_test,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
        disp_groups = c('case' = TRUE, 'control' = TRUE),
        local_window_size = 0,
        t_approx = TRUE,
        n_cores = 1)

    expect_equal(2+2, 4)
})

test_that('Check dropped loci message', {
    expect_message(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = FALSE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'loci were dropped due to insufficient degrees',
        fixed = TRUE
    )
})

test_that('Test 2', {
    diff_gr = diff_methylsig(
        bs = small_test,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
        disp_groups = c('case' = TRUE, 'control' = FALSE),
        local_window_size = 0,
        t_approx = TRUE,
        n_cores = 1)

    expect_equal(2+2, 4)
})

test_that('Test 3', {
    diff_gr = diff_methylsig(
        bs = small_test,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
        disp_groups = c('case' = FALSE, 'control' = TRUE),
        local_window_size = 0,
        t_approx = TRUE,
        n_cores = 1)

    expect_equal(2+2, 4)
})

test_that('Test 4', {
    diff_gr = diff_methylsig(
        bs = small_test_tile,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
        disp_groups = c('case' = FALSE, 'control' = TRUE),
        local_window_size = 0,
        t_approx = TRUE,
        n_cores = 1)

    expect_equal(2+2, 4)
})

test_that('Test 5', {
    diff_gr = diff_methylsig(
        bs = small_test,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
        disp_groups = c('case' = FALSE, 'control' = TRUE),
        local_window_size = 0,
        t_approx = FALSE,
        n_cores = 1)

    expect_equal(2+2, 4)
})

test_that('Test 6', {
    diff_gr = diff_methylsig(
        bs = small_test,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
        disp_groups = c('case' = FALSE, 'control' = TRUE),
        local_window_size = 50,
        t_approx = TRUE,
        n_cores = 1)

    expect_equal(2+2, 4)
})
