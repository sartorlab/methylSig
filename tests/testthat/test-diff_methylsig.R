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

    expect_true(is(diff_gr, 'GRanges'))
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

    expect_true(is(diff_gr, 'GRanges'))
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

    expect_true(is(diff_gr, 'GRanges'))
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

    expect_true(is(diff_gr, 'GRanges'))
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

    expect_true(is(diff_gr, 'GRanges'))
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

    expect_true(is(diff_gr, 'GRanges'))
})

test_that('Untested loci are dropped', {
    # Tiles from position 1 onward, most with no coverage, and so untested
    diff_gr = suppressMessages(diff_methylsig(
        bs = tile_by_windows(bs = small_test, win_size = 100),
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
        disp_groups = c('case' = TRUE, 'control' = TRUE),
        local_window_size = 0,
        t_approx = TRUE,
        n_cores = 1))

    expect_gt(length(diff_gr), 0)
    expect_false(anyNA(diff_gr$pvalue))
    expect_true(all(diff_gr$df > 3))
})

test_that('One sample in a group', {
    one_case = small_test[, c(1, 4, 5, 6)]
    expect_equal(sum(bsseq::pData(one_case)$Type == 'cancer'), 1)

    diff_gr = suppressMessages(diff_methylsig(
        bs = one_case,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
        disp_groups = c('case' = FALSE, 'control' = TRUE),
        local_window_size = 0,
        t_approx = TRUE,
        n_cores = 1))

    expect_true(is(diff_gr, 'GRanges'))
    expect_gt(length(diff_gr), 0)
})

test_that('Too few samples to estimate dispersion', {
    expect_error(
        diff_methylsig(
            bs = small_test[, c(1, 4, 5, 6)],
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = FALSE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'Too few samples to estimate dispersion: disp_groups has 1 samples, and at least 3 are needed.',
        fixed = TRUE
    )
    expect_error(
        diff_methylsig(
            bs = small_test[, c(1, 2, 4)],
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 0,
            t_approx = TRUE,
            n_cores = 1),
        'Too few samples to estimate dispersion: disp_groups has 3 samples, and at least 4 are needed.',
        fixed = TRUE
    )
})

test_that('Local information stays within a chromosome', {
    # Copy small_test to chr22, starting 10bp after the last chr21 locus, so
    # the first chr22 loci are within local_window_size of the last chr21
    # loci by position alone
    gr21 = GenomicRanges::granges(small_test)
    shift = max(BiocGenerics::start(gr21)) - min(BiocGenerics::start(gr21)) + 10
    gr22 = GenomicRanges::GRanges('chr22', IRanges::IRanges(BiocGenerics::start(gr21) + shift, width = 1))
    two_chrom = bsseq::BSseq(
        gr = c(gr21, gr22),
        M = rbind(as.matrix(bsseq::getCoverage(small_test, type = 'M')), as.matrix(bsseq::getCoverage(small_test, type = 'M'))),
        Cov = rbind(as.matrix(bsseq::getCoverage(small_test, type = 'Cov')), as.matrix(bsseq::getCoverage(small_test, type = 'Cov'))),
        pData = bsseq::pData(small_test),
        sampleNames = colnames(small_test))

    run = function(bs) {
        suppressMessages(diff_methylsig(
            bs = bs,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 200,
            t_approx = TRUE,
            n_cores = 1))
    }
    alone = run(small_test)
    both = run(two_chrom)
    both21 = both[GenomicRanges::seqnames(both) == 'chr21']

    # fdr depends on the number of tests, so compare everything else
    cols = setdiff(names(S4Vectors::mcols(alone)), 'fdr')
    expect_equal(BiocGenerics::start(both21), BiocGenerics::start(alone))
    expect_equal(as.data.frame(S4Vectors::mcols(both21)[, cols]), as.data.frame(S4Vectors::mcols(alone)[, cols]))
})

test_that('local_disp and local_meth type checks', {
    for (arg in c('local_disp', 'local_meth')) {
        args = list(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 50)
        args[[arg]] = NA
        expect_error(do.call(diff_methylsig, args), sprintf('%s must be TRUE/FALSE.', arg), fixed = TRUE)
    }
    expect_error(
        diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = 50,
            local_disp = FALSE,
            local_meth = FALSE),
        'local_window_size > 0 uses local information, so local_disp or local_meth must be TRUE.',
        fixed = TRUE
    )
})

test_that('local_disp and local_meth use local information separately', {
    run = function(local_window_size, local_disp = TRUE, local_meth = TRUE) {
        suppressMessages(diff_methylsig(
            bs = small_test,
            group_column = 'Type',
            comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
            disp_groups = c('case' = TRUE, 'control' = TRUE),
            local_window_size = local_window_size,
            t_approx = TRUE,
            n_cores = 1,
            local_disp = local_disp,
            local_meth = local_meth))
    }
    none = run(0)
    both = run(50)
    disp_only = run(50, local_meth = FALSE)
    meth_only = run(50, local_disp = FALSE)

    # Dispersion and degrees of freedom come from the dispersion data
    expect_equal(disp_only$disp_est, both$disp_est)
    expect_equal(disp_only$df, both$df)
    expect_equal(meth_only$disp_est, none$disp_est)
    expect_equal(meth_only$df, none$df)

    # Group methylation and the statistic come from the methylation data
    expect_false(isTRUE(all.equal(disp_only$log_lik_ratio, both$log_lik_ratio)))
    expect_false(isTRUE(all.equal(meth_only$log_lik_ratio, none$log_lik_ratio)))
})

test_that('Local log likelihood is the weighted sum over loci', {
    local_c = matrix(c(3, 0, 7, 5, 2, 9), nrow = 3)
    local_t = matrix(c(4, 6, 0, 1, 8, 2), nrow = 3)
    weights = c(0.5, 1, 0.25)
    mu = 0.4
    phi = 12

    per_locus = vapply(seq(nrow(local_c)), function(i) {
        .log_likelihood(
            mu = mu,
            phi = phi,
            local_c = local_c[i, , drop = FALSE],
            local_t = local_t[i, , drop = FALSE],
            weight = 1)
    }, 1)

    expect_equal(
        .log_likelihood(mu, phi, local_c, local_t, weight = rep(1, 3)),
        sum(per_locus))
    expect_equal(
        .log_likelihood(mu, phi, local_c, local_t, weight = weights),
        sum(weights * per_locus))
})

test_that('Too few samples to estimate dispersion without local dispersion', {
    three = small_test[, c(1, 2, 4)]
    args = list(
        bs = three,
        group_column = 'Type',
        comparison_groups = c('case' = 'cancer', 'control' = 'normal'),
        disp_groups = c('case' = TRUE, 'control' = TRUE),
        local_window_size = 50,
        t_approx = TRUE,
        n_cores = 1)

    # Local dispersion can add degrees of freedom, so 3 samples are enough
    expect_true(is(suppressMessages(do.call(diff_methylsig, args)), 'GRanges'))

    args$local_disp = FALSE
    expect_error(
        do.call(diff_methylsig, args),
        'Too few samples to estimate dispersion: disp_groups has 3 samples, and at least 4 are needed.',
        fixed = TRUE
    )
})
