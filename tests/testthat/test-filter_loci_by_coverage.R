test_that('BSseq class check', {
    expect_error(
        filter_loci_by_coverage(5),
        'bs must be class BSseq'
    )
})

test_that('min_count numeric check', {
    expect_error(
        filter_loci_by_coverage(bsseq_stranded, min_count = 'a'),
        'min_count must be an integer'
    )
})

test_that('max_count numeric check', {
    expect_error(
        filter_loci_by_coverage(bsseq_stranded, max_count = 'a'),
        'max_count must be an integer'
    )
})

test_that('min_count less than max_count check', {
    expect_error(
        filter_loci_by_coverage(bsseq_stranded, min_count = 600 ),
        'min_count not less than max_count'
    )
})

test_that('correct set to 0 check', {
    bs = filter_loci_by_coverage(bsseq_stranded, min_count = 10, max_count = 500)
    bs_cov = bsseq::getCoverage(bs, type = 'Cov')
    bs_meth = bsseq::getCoverage(bs, type = 'M')

    expect_equivalent(bs_cov[1,'test1'], 0)
    expect_equivalent(bs_cov[2,'test1'], 0)
    expect_equivalent(bs_cov[7,'test2'], 0)
    expect_equivalent(bs_cov[8,'test2'], 0)
    expect_equivalent(bs_cov[10,'test1'], 0)
    expect_equivalent(bs_cov[11,'test1'], 0)

    expect_equivalent(bs_meth[1,'test1'], 0)
    expect_equivalent(bs_meth[2,'test1'], 0)
    expect_equivalent(bs_meth[7,'test2'], 0)
    expect_equivalent(bs_meth[8,'test2'], 0)
    expect_equivalent(bs_meth[10,'test1'], 0)
    expect_equivalent(bs_meth[11,'test1'], 0)
})
