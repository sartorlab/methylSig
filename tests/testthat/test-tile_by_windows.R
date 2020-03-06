test_that('bs missing check', {
    expect_error(
        tile_by_windows(),
        'Must pass bs as a BSseq object',
        fixed = TRUE
    )
})

test_that('bs class check', {
    expect_error(
        tile_by_windows(bs = '5', win_size = 200),
        'bs must be class BSseq',
        fixed = TRUE
    )
})

test_that('win_size class check', {
    expect_error(
        tile_by_windows(bs = bsseq_stranded, win_size = TRUE),
        'win_size must be an integer',
        fixed = TRUE
    )
})

test_that('correct tiling stranded win25', {
    test = tile_by_windows(bs = bsseq_stranded, win_size = 25)

    # NOTE, it is not sufficient to do
    # expect_equivalent(test, bsseq_stranded_tiled1) because the test
    # expect_equivalent(test, bsseq_destranded_tiled1) does not throw an
    # error but testing equivalance of the Cov or M matrices of test and
    # bsseq_destranded_tiled1 will. Consequently, testing at this
    # level is necessary

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_stranded_win25, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_stranded_win25, type = 'M')
    )

    expect_true(
        all(is.na(seqlengths(test)))
    )
})

test_that('correct tiling destranded win25', {
    test = tile_by_windows(bs = bsseq_destranded, win_size = 25)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_destranded_win25, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_destranded_win25, type = 'M')
    )

    expect_true(
        all(is.na(seqlengths(test)))
    )
})

test_that('correct tiling multichrom multichrom25', {
    test = tile_by_windows(bs = bsseq_multichrom, win_size = 25)

    expect_equivalent(
        bsseq::getCoverage(test, type = 'Cov'),
        bsseq::getCoverage(bsseq_multichrom_win25, type = 'Cov')
    )

    expect_equivalent(
        bsseq::getCoverage(test, type = 'M'),
        bsseq::getCoverage(bsseq_multichrom_win25, type = 'M')
    )

    expect_true(
        all(is.na(seqlengths(test)))
    )
})
