context('Test methylSigTile')

################################################################################

# Tiles on the same chromosomes as data
tiles_df_samechr = read.table(system.file('extdata','test_tiles.txt', package='methylSig'), header = T, sep = '\t', as.is = T)
tiles_gr_samechr_noseqinfo = GenomicRanges::makeGRangesFromDataFrame(tiles_df_samechr)
tiles_gr_samechr_fullseqinfo = GenomicRanges::makeGRangesFromDataFrame(tiles_df_samechr, seqinfo = GenomeInfoDb::Seqinfo(genome='hg19'))

# Tiles with some tiles on different chromosomes as data
tiles_df_extrachr = read.table(system.file('extdata','test_tiles2.txt', package='methylSig'), header = T, sep = '\t', as.is = T)
tiles_gr_extrachr_noseqinfo = GenomicRanges::makeGRangesFromDataFrame(tiles_df_extrachr)
tiles_gr_extrachr_fullseqinfo = GenomicRanges::makeGRangesFromDataFrame(tiles_df_extrachr, seqinfo = GenomeInfoDb::Seqinfo(genome='hg19'))

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

################################################################################

truth1_meth = matrix(c(34,10,0,300,87,434), nrow = 3, ncol = 2, byrow = TRUE)
truth1_cov = matrix(c(34,10,0,300,96,458), nrow = 3, ncol = 2, byrow = TRUE)

test_that('Test tileGenome tiling', {
    tiled_data = methylSigTile(
        meth = data,
        tiles = NULL,
        win.size = 200)

    expect_true(all(truth1_meth == as.matrix(bsseq::getCoverage(tiled_data, type='M'))))
    expect_true(all(truth1_cov == as.matrix(bsseq::getCoverage(tiled_data, type='Cov'))))
    expect_true(unique(GenomeInfoDb::genome(tiled_data)) == 'hg19')
})

truth2_meth = matrix(c(34,10,87,734), nrow = 2, ncol = 2, byrow = TRUE)
truth2_cov = matrix(c(34,10,96,758), nrow = 2, ncol = 2, byrow = TRUE)

test_that('Test data.frame tiling with matching chromosomes', {
    tiled_data = methylSigTile(
        meth = data,
        tiles = tiles_df_samechr,
        win.size = 200)

    expect_true(all(truth2_meth == as.matrix(bsseq::getCoverage(tiled_data, type='M'))))
    expect_true(all(truth2_cov == as.matrix(bsseq::getCoverage(tiled_data, type='Cov'))))
    expect_true(unique(GenomeInfoDb::genome(tiled_data)) == 'hg19')
})

test_that('Test GRanges tiling with matching chromosomes and no seqinfo', {
    tiled_data = methylSigTile(
        meth = data,
        tiles = tiles_gr_samechr_noseqinfo,
        win.size = 200)

    expect_true(all(truth2_meth == as.matrix(bsseq::getCoverage(tiled_data, type='M'))))
    expect_true(all(truth2_cov == as.matrix(bsseq::getCoverage(tiled_data, type='Cov'))))
    expect_true(unique(GenomeInfoDb::genome(tiled_data)) == 'hg19')
})

test_that('Test GRanges tiling with matching chromosomes and full seqinfo', {
    tiled_data = methylSigTile(
        meth = data,
        tiles = tiles_gr_samechr_fullseqinfo,
        win.size = 200)

    expect_true(all(truth2_meth == as.matrix(bsseq::getCoverage(tiled_data, type='M'))))
    expect_true(all(truth2_cov == as.matrix(bsseq::getCoverage(tiled_data, type='Cov'))))
    expect_true(unique(GenomeInfoDb::genome(tiled_data)) == 'hg19')
})

test_that('Test data.frame tiling with extra chromosomes', {
    tiled_data = methylSigTile(
        meth = data,
        tiles = tiles_df_extrachr,
        win.size = 200)

    expect_true(all(truth2_meth == as.matrix(bsseq::getCoverage(tiled_data, type='M'))))
    expect_true(all(truth2_cov == as.matrix(bsseq::getCoverage(tiled_data, type='Cov'))))
    expect_true(unique(GenomeInfoDb::genome(tiled_data)) == 'hg19')
})

test_that('Test GRanges tiling with extra chromosomes and no seqinfo', {
    tiled_data = methylSigTile(
        meth = data,
        tiles = tiles_gr_extrachr_noseqinfo,
        win.size = 200)

    expect_true(all(truth2_meth == as.matrix(bsseq::getCoverage(tiled_data, type='M'))))
    expect_true(all(truth2_cov == as.matrix(bsseq::getCoverage(tiled_data, type='Cov'))))
    expect_true(unique(GenomeInfoDb::genome(tiled_data)) == 'hg19')
})

test_that('Test GRanges tiling with extra chromosomes and full seqinfo', {
    tiled_data = methylSigTile(
        meth = data,
        tiles = tiles_gr_extrachr_fullseqinfo,
        win.size = 200)

    expect_true(all(truth2_meth == as.matrix(bsseq::getCoverage(tiled_data, type='M'))))
    expect_true(all(truth2_cov == as.matrix(bsseq::getCoverage(tiled_data, type='Cov'))))
    expect_true(unique(GenomeInfoDb::genome(tiled_data)) == 'hg19')
})

test_that('Test GRanges tiling warning', {
    expect_message(
        methylSigTile(
            meth = data,
            tiles = tiles_gr_samechr_noseqinfo,
            win.size = 200),
        'genome of the GRanges')
})
