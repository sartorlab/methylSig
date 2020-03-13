test_that('bs missing check', {
    expect_error(
        filter_loci_by_location(),
        'Must pass bs as a BSseq object',
        fixed = TRUE
    )
})

test_that('gr missing check', {
    expect_error(
        filter_loci_by_location(bs = bsseq_stranded),
        'Must pass gr as a GRanges object',
        fixed = TRUE
    )
})

test_that('bs class check', {
    expect_error(
        filter_loci_by_location(bs = '5', gr = gr_tiles1),
        'bs must be class BSseq',
        fixed = TRUE
    )
})

test_that('gr class check', {
    expect_error(
        filter_loci_by_location(bs = bsseq_stranded, gr = '5'),
        'gr must be class GRanges',
        fixed = TRUE
    )
})

#####################################

test_that('correct filtering gr1', {
    test = filter_loci_by_location(bs = bsseq_stranded, gr = gr_tiles1)

    expect_equivalent(
        as.matrix(bsseq::getCoverage(test, type = 'Cov')),
        as.matrix(bsseq::getCoverage(filter_loc_tiles1, type = 'Cov'))
    )

    expect_equivalent(
        as.matrix(bsseq::getCoverage(test, type = 'M')),
        as.matrix(bsseq::getCoverage(filter_loc_tiles1, type = 'M'))
    )
})

test_that('correct filtering gr3', {
    expect_error(
        filter_loci_by_location(bs = bsseq_stranded, gr = gr_tiles3),
        'All loci in bs were removed by gr',
        fixed = TRUE
    )
})

test_that('correct filtering gr4', {
    test = filter_loci_by_location(bs = bsseq_stranded, gr = gr_tiles4)

    expect_equivalent(
        as.matrix(bsseq::getCoverage(test, type = 'Cov')),
        as.matrix(bsseq::getCoverage(filter_loc_tiles4, type = 'Cov'))
    )

    expect_equivalent(
        as.matrix(bsseq::getCoverage(test, type = 'M')),
        as.matrix(bsseq::getCoverage(filter_loc_tiles4, type = 'M'))
    )
})

test_that('correct filtering gr5', {
    test = filter_loci_by_location(bs = bsseq_stranded, gr = gr_tiles5)

    expect_equivalent(
        as.matrix(bsseq::getCoverage(test, type = 'Cov')),
        as.matrix(bsseq::getCoverage(filter_loc_tiles5, type = 'Cov'))
    )

    expect_equivalent(
        as.matrix(bsseq::getCoverage(test, type = 'M')),
        as.matrix(bsseq::getCoverage(filter_loc_tiles5, type = 'M'))
    )
})
