test_that('bs missing check', {
    expect_error(
        tile_by_regions(),
        'Must pass bs as a BSseq object'
    )
})

test_that('gr missing check', {
    expect_error(
        tile_by_regions(bs = bsseq_cov_with_strand),
        'Must pass gr as a GRanges object'
    )
})

test_that('bs class check', {
    expect_error(
        tile_by_regions(bs = '5', gr = tile_regions_gr1),
        'bs must be class BSseq'
    )
})

test_that('gr class check', {
    expect_error(
        tile_by_regions(bs = bsseq_cov_with_strand, gr = '5'),
        'gr must be class GRanges'
    )
})


test_that('correct tiling with strand gr1', {
    test = tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr1)

    # NOTE, it is not sufficient to do
    # expect_equivalent(test, bsseq_with_strand_tile_gr1) because the test
    # expect_equivalent(test, bsseq_without_strand_tile_gr1) does not throw an
    # error but testing equivalance of the Cov or M matrices of test and
    # bsseq_without_strand_tile_gr1 will. Consequently, testing at this
    # level is necessary

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_with_strand_tile_gr1, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_with_strand_tile_gr1, type = 'M')
    )
})

test_that('correct tiling with strand gr2', {
    test = tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr2)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_with_strand_tile_gr2, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_with_strand_tile_gr2, type = 'M')
    )
})

test_that('correct tiling with strand gr3', {
    test = tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr3)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_with_strand_tile_gr3, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_with_strand_tile_gr3, type = 'M')
    )
})

test_that('correct tiling with strand gr4', {
    test = tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr4)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_with_strand_tile_gr4, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_with_strand_tile_gr4, type = 'M')
    )
})

test_that('error tiling with strand gr5', {
    expect_error(
        tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr5),
        'No regions overlap between bs and gr'
    )
})

test_that('correct tiling without strand gr1', {
    test = tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr1)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_without_strand_tile_gr1, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_without_strand_tile_gr1, type = 'M')
    )
})

test_that('correct tiling without strand gr2', {
    test = tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr2)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_without_strand_tile_gr2, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_without_strand_tile_gr2, type = 'M')
    )
})

test_that('correct tiling without strand gr3', {
    test = tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr3)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_without_strand_tile_gr3, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_without_strand_tile_gr3, type = 'M')
    )
})

test_that('correct tiling without strand gr4', {
    test = tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr4)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_without_strand_tile_gr4, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_without_strand_tile_gr4, type = 'M')
    )
})

test_that('error tiling without strand gr5', {
    expect_error(
        tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr5),
        'No regions overlap between bs and gr'
    )
})
