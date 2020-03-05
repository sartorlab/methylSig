data(BS.cancer.ex, package = 'bsseqData')

small_test = BS.cancer.ex[1:10]

expected_cov_cancer2_normal2 = bsseq::getCoverage(small_test, type = 'Cov')[c(3,4,7,8,9,10), ]
expected_meth_cancer2_normal2 = bsseq::getCoverage(small_test, type = 'M')[c(3,4,7,8,9,10), ]

expected_cov_cancer2_normal3 = bsseq::getCoverage(small_test, type = 'Cov')[c(3,8,9), ]
expected_meth_cancer2_normal3 = bsseq::getCoverage(small_test, type = 'M')[c(3,8,9), ]

expected_cov_cancer3_normal3 = bsseq::getCoverage(small_test, type = 'Cov')[c(3,8), ]
expected_meth_cancer3_normal3 = bsseq::getCoverage(small_test, type = 'M')[c(3,8), ]

test_that('bs missing check', {
    expect_error(
        filter_loci_by_group_coverage(),
        'Must pass bs as a BSseq object'
    )
})

test_that('group_column missing check', {
    expect_error(
        filter_loci_by_group_coverage(bs = small_test),
        'Must pass group_column as a character string'
    )
})

test_that('min_samples_per_group missing check', {
    expect_error(
        filter_loci_by_group_coverage(bs = small_test, group_column = 'Type'),
        'Must pass min_samples_per_group as a named integer vector'
    )
})

test_that('bs type check', {
    expect_error(
        filter_loci_by_group_coverage(
            bs = 'a',
            group_column = 'Type',
            c('cancer' = 3, 'normal' = 3)),
        'bs must be class BSseq'
    )
})

test_that('group_column type check', {
    expect_error(
        filter_loci_by_group_coverage(
            bs = small_test,
            group_column = 6,
            c('cancer' = 3, 'normal' = 3)),
        'group_column must be a character string'
    )
})

test_that('min_samples_per_group type check', {
    expect_error(
        filter_loci_by_group_coverage(
            bs = small_test,
            group_column = 'Type',
            c('cancer' = 'a', 'normal' = 3)),
        'min_samples_per_group must be a named integer vector'
    )
})

test_that('Valid group_column name check', {
    expect_error(
        filter_loci_by_group_coverage(
            bs = small_test,
            group_column = 'blue',
            c('cancer' = 3, 'normal' = 3)),
        'group_column: blue not in column names of pData' # () seem to be a problem
    )
})

test_that('Valid factor name check', {
    expect_error(
        filter_loci_by_group_coverage(
            bs = small_test,
            group_column = 'Type',
            c('blue' = 3, 'normal' = 3)),
        'Not all names of min_sample_per_group are in group_column'
    )
})

test_that('All loci removed check', {
    expect_error(
        filter_loci_by_group_coverage(
            bs = small_test,
            group_column = 'Type',
            c('cancer' = 4, 'normal' = 4)),
        'Thresholds for the following groups were too strict'
    )
})

test_that('Test cancer 2 normal 2', {
    test = filter_loci_by_group_coverage(
        bs = small_test,
        group_column = 'Type',
        c('cancer' = 2, 'normal' = 2))

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        expected_cov_cancer2_normal2
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        expected_meth_cancer2_normal2
    )
})

test_that('Test cancer 2 normal 3', {
    test = filter_loci_by_group_coverage(
        bs = small_test,
        group_column = 'Type',
        c('cancer' = 2, 'normal' = 3))

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        expected_cov_cancer2_normal3
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        expected_meth_cancer2_normal3
    )
})

test_that('Test cancer 3 normal 3', {
    test = filter_loci_by_group_coverage(
        bs = small_test,
        group_column = 'Type',
        c('cancer' = 3, 'normal' = 3))

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        expected_cov_cancer3_normal3
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        expected_meth_cancer3_normal3
    )
})
