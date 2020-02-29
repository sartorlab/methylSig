test_that('bs missing check', {
    expect_error(
        tile_by_regions(),
        'Must pass bs as a BSseq object'
    )
})

test_that('gr missing check', {
    expect_error(
        tile_by_regions(bs = bsseq_stranded),
        'Must pass gr as a GRanges object'
    )
})

test_that('bs class check', {
    expect_error(
        tile_by_regions(bs = '5', gr = gr_tiles1),
        'bs must be class BSseq'
    )
})

test_that('gr class check', {
    expect_error(
        tile_by_regions(bs = bsseq_stranded, gr = '5'),
        'gr must be class GRanges'
    )
})


test_that('correct tiling stranded gr1', {
    test = tile_by_regions(bs = bsseq_stranded, gr = gr_tiles1)

    # NOTE, it is not sufficient to do
    # expect_equivalent(test, bsseq_stranded_tiled1) because the test
    # expect_equivalent(test, bsseq_destranded_tiled1) does not throw an
    # error but testing equivalance of the Cov or M matrices of test and
    # bsseq_destranded_tiled1 will. Consequently, testing at this
    # level is necessary

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_stranded_tiled1, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_stranded_tiled1, type = 'M')
    )
})

test_that('correct tiling stranded gr2', {
    test = tile_by_regions(bs = bsseq_stranded, gr = gr_tiles2)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_stranded_tiled2, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_stranded_tiled2, type = 'M')
    )
})

test_that('correct tiling stranded gr3', {
    test = tile_by_regions(bs = bsseq_stranded, gr = gr_tiles3)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_stranded_tiled3, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_stranded_tiled3, type = 'M')
    )
})

test_that('correct tiling stranded gr4', {
    test = tile_by_regions(bs = bsseq_stranded, gr = gr_tiles4)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_stranded_tiled4, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_stranded_tiled4, type = 'M')
    )
})

test_that('error tiling stranded gr5', {
    expect_error(
        tile_by_regions(bs = bsseq_stranded, gr = gr_tiles5),
        'No regions overlap between bs and gr'
    )
})

test_that('correct tiling destranded gr1', {
    test = tile_by_regions(bs = bsseq_destranded, gr = gr_tiles1)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_destranded_tiled1, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_destranded_tiled1, type = 'M')
    )
})

test_that('correct tiling destranded gr2', {
    test = tile_by_regions(bs = bsseq_destranded, gr = gr_tiles2)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_destranded_tiled2, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_destranded_tiled2, type = 'M')
    )
})

test_that('correct tiling destranded gr3', {
    test = tile_by_regions(bs = bsseq_destranded, gr = gr_tiles3)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_destranded_tiled3, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_destranded_tiled3, type = 'M')
    )
})

test_that('correct tiling destranded gr4', {
    test = tile_by_regions(bs = bsseq_destranded, gr = gr_tiles4)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_destranded_tiled4, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_destranded_tiled4, type = 'M')
    )
})

test_that('error tiling destranded gr5', {
    expect_error(
        tile_by_regions(bs = bsseq_destranded, gr = gr_tiles5),
        'No regions overlap between bs and gr'
    )
})
