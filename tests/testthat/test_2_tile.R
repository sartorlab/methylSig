context('Test methylSigTile')

################################################################################

tiles_df = read.table(system.file('extdata','test_tiles.txt', package='methylSig'), header = T, sep = '\t', as.is = T)

tiles_gr = makeGRangesFromDataFrame(tiles_df)

files = c(system.file('extdata', 'test_1.txt', package='methylSig'),
    system.file('extdata', 'test_2.txt', package='methylSig'))

sample_names = gsub('.txt', '', basename(files))

pData = data.frame(
    Sample_Names = sample_names,
    Group = relevel(factor(c(1,0)), ref = '0'),
    Note = c("Hello", "Goodbye"),
    row.names = sample_names,
    stringsAsFactors = FALSE)

data = methylSigReadData(
    fileList = files,
    pData = pData,
    assembly = 'hg19',
    destranded = TRUE,
    maxCount = 500,
    minCount = 10,
    filterSNPs = TRUE,
    num.cores = 1,
    fileType = 'cytosineReport')

tiles_gr_seqinfo = tiles_gr
seqinfo(tiles_gr_seqinfo) = seqinfo(data)

################################################################################

truth1_meth = matrix(c(34,10,0,300,87,434), nrow = 3, ncol = 2, byrow = TRUE)
truth1_cov = matrix(c(34,10,0,300,96,458), nrow = 3, ncol = 2, byrow = TRUE)

test_that('Test tileGenome tiling', {
    tiled_data1 = methylSigTile(
        meth = data,
        tiles = NULL,
        win.size = 200)

    expect_true(all(truth1_meth == as.matrix(bsseq::getCoverage(tiled_data1, type='M'))))
    expect_true(all(truth1_cov == as.matrix(bsseq::getCoverage(tiled_data1, type='Cov'))))
})

truth2_meth = matrix(c(34,10,87,734), nrow = 2, ncol = 2, byrow = TRUE)
truth2_cov = matrix(c(34,10,96,758), nrow = 2, ncol = 2, byrow = TRUE)

test_that('Test data.frame tiling', {
    tiled_data2 = methylSigTile(
        meth = data,
        tiles = tiles_df,
        win.size = 200)

    expect_true(all(truth2_meth == as.matrix(bsseq::getCoverage(tiled_data2, type='M'))))
    expect_true(all(truth2_cov == as.matrix(bsseq::getCoverage(tiled_data2, type='Cov'))))
})

test_that('Test GRanges tiling warning', {
    expect_warning(
        methylSigTile(
            meth = data,
            tiles = tiles_gr,
            win.size = 200),
        'genome of the GRanges')
})

test_that('Test GRanges tiling', {
    tiled_data3 = suppressWarnings(methylSigTile(
        meth = data,
        tiles = tiles_gr,
        win.size = 200))

    expect_true(all(truth2_meth == as.matrix(bsseq::getCoverage(tiled_data3, type='M'))))
    expect_true(all(truth2_cov == as.matrix(bsseq::getCoverage(tiled_data3, type='Cov'))))
})

test_that('Test GRanges with seqinfo tiling', {
    tiled_data4 = methylSigTile(
        meth = data,
        tiles = tiles_gr_seqinfo,
        win.size = 200)

    expect_true(all(truth2_meth == as.matrix(bsseq::getCoverage(tiled_data4, type='M'))))
    expect_true(all(truth2_cov == as.matrix(bsseq::getCoverage(tiled_data4, type='Cov'))))
})
