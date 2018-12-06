context('Test methylSigCalc')

utils::data(sample_data, package = 'methylSig')

test_that('Test CpG no local both dispersion', {
    result_calc = methylSigCalc(
        meth = meth,
        comparison = 'DR_vs_DS',
        dispersion = 'both',
        local.info = FALSE,
        local.winsize = 200,
        min.per.group = c(3,3),
        weightFunc = methylSig_weightFunc,
        T.approx = TRUE,
        num.cores = 1)

    expect_true(is(result_calc, 'GRanges'))
    expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
})

# test_that('Test CpG no local DR dispersion', {
#     result_calc = methylSigCalc(
#         meth = meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'DR',
#         local.info = FALSE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })

# test_that('Test CpG no local DS dispersion', {
#     result_calc = methylSigCalc(
#         meth = meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'DS',
#         local.info = FALSE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })

# test_that('Test CpG local both dispersion', {
#     result_calc = methylSigCalc(
#         meth = meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'both',
#         local.info = TRUE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })

# test_that('Test CpG local DR dispersion', {
#     result_calc = methylSigCalc(
#         meth = meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'DR',
#         local.info = TRUE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })

# test_that('Test CpG local DS dispersion', {
#     result_calc = methylSigCalc(
#         meth = meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'DS',
#         local.info = TRUE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })

test_that('Test tiled no local both dispersion', {
    result_calc = methylSigCalc(
        meth = tiled_meth,
        comparison = 'DR_vs_DS',
        dispersion = 'both',
        local.info = FALSE,
        local.winsize = 200,
        min.per.group = c(3,3),
        weightFunc = methylSig_weightFunc,
        T.approx = TRUE,
        num.cores = 1)

    expect_true(is(result_calc, 'GRanges'))
    expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
})

# test_that('Test tiled no local DR dispersion', {
#     result_calc = methylSigCalc(
#         meth = tiled_meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'DR',
#         local.info = FALSE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })

# test_that('Test tiled no local DS dispersion', {
#     result_calc = methylSigCalc(
#         meth = tiled_meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'DS',
#         local.info = FALSE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })

# test_that('Test tiled local both dispersion', {
#     result_calc = methylSigCalc(
#         meth = tiled_meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'both',
#         local.info = TRUE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })

# test_that('Test tiled local DR dispersion', {
#     result_calc = methylSigCalc(
#         meth = tiled_meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'DR',
#         local.info = TRUE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })

# test_that('Test tiled local DS dispersion', {
#     result_calc = methylSigCalc(
#         meth = tiled_meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'DS',
#         local.info = TRUE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })


# test_that('Test specified tiles', {
#     cpg_islands = annotatr::build_annotations(genome = 'hg19', annotations = 'hg19_cpg_islands')
#     cpg_island_meth = methylSigTile(meth, tiles = cpg_islands)
#
#     result_calc = methylSigCalc(
#         meth = cpg_island_meth,
#         comparison = 'DR_vs_DS',
#         dispersion = 'both',
#         local.info = FALSE,
#         local.winsize = 200,
#         min.per.group = c(3,3),
#         weightFunc = methylSig_weightFunc,
#         T.approx = TRUE,
#         num.cores = 1)
#
#     expect_true(is(result_calc, 'GRanges'))
#     expect_match(S4Vectors::metadata(result_calc)$method, 'methylSigCalc')
# })
