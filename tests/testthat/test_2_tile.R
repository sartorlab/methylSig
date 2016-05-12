context('Test methylSigTile')

################################################################################
# Test with error message

files = c(
  system.file('extdata', 'test_1.txt', package='methylSig'),
  system.file('extdata', 'test_2.txt', package='methylSig'))

data_error_tiled = methylSigReadData(
  fileList = files,
  sample.ids = c('test_1','test_2'),
  assembly = 'hg19',
  pipeline = 'bismark and methylKit',
  header = TRUE,
  context = 'CpG',
  resolution = "region",
  treatment = c(1,0),
  destranded = TRUE,
  maxCount = 500,
  minCount = 10,
  filterSNPs = FALSE,
  num.cores = 1,
  quiet = TRUE)

test_that('Already tiled data throws error', {
  expect_error(methylSigTile(data_error_tiled), 'Object has already been tiled')
})

################################################################################
# Test window tiling with default settings (except header)

data_default = methylSigReadData(
  fileList = files,
  sample.ids = c('test_1','test_2'),
  assembly = 'hg19',
  pipeline = 'bismark and methylKit',
  header = FALSE,
  context = 'CpG',
  resolution = "base",
  treatment = c(1,0),
  destranded = TRUE,
  maxCount = 500,
  minCount = 10,
  filterSNPs = FALSE,
  num.cores = 1,
  quiet = TRUE)

tiled_windows = methylSigTile(meth = data_default, win.size = 25)

test_that('Test tiling by window produces 25bp windows except for 1', {
  expect_equal( length(which(names(table(tiled_windows@data.end - tiled_windows@data.start + 1)) != '25')), expected = 1)
})

test_that('Test tiling does not exceed chromosome bound',{
  # While methylSig is genome ignorant, we can avoid exceeding chromosome
  # bounds by extending tiles only to max(meth@data.end)
  expect_equal( tiled_windows@data.end[length(tiled_windows@data.end)] > 48129895, expected = FALSE)
})

test_that('Test that tiling aggregates properly', {
  expect_equal( tiled_windows@data.coverage[which(tiled_windows@data.start == 43053285), 1], expected = 96)
  expect_equal( tiled_windows@data.coverage[which(tiled_windows@data.start == 43053285), 2], expected = 274)
})

################################################################################
# Test annotation tiling with default parameters

cgi_file = system.file('extdata', 'test_annotation.txt', package='methylSig')
cgis = read.table(cgi_file, header=F, sep='\t', stringsAsFactors=F)
colnames(cgis) = c('chr','start','end','name')

tiled_cgis = methylSigTile(meth = data_default, tiles = cgis)

test_that('Test tiling by test_annotation aggregates properly',{
  expect_equal( tiled_cgis@data.coverage[which(tiled_cgis@data.start == 43045070), 1], expected = 44)
  expect_equal( tiled_cgis@data.coverage[which(tiled_cgis@data.start == 43045070), 2], expected = 0)
  expect_equal( tiled_cgis@data.coverage[which(tiled_cgis@data.start == 43052350), 1], expected = 116)
  expect_equal( tiled_cgis@data.coverage[which(tiled_cgis@data.start == 43052350), 2], expected = 207)
  expect_equal( tiled_cgis@data.coverage[which(tiled_cgis@data.start == 43053316), 1], expected = 0)
  expect_equal( tiled_cgis@data.coverage[which(tiled_cgis@data.start == 43053316), 2], expected = 176)
})
