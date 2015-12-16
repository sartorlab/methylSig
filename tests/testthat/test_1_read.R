context('Test methylSigReadData')

################################################################################
# Test with default settings (except header)

files = c(
  system.file('extdata', 'test_1.txt', package='methylSig'),
  system.file('extdata', 'test_2.txt', package='methylSig'))

data_defaults = methylSigReadData(
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

test_that('Test minCount removals after destrand',{
  expect_equal( 43053323 %in% data_defaults@data.start, expected = FALSE )
})

test_that('Test maxCount removals after destrand',{
  expect_equal( 43052356 %in% data_defaults@data.start, expected = FALSE )
})

test_that('Test invalid frequency removals',{
  expect_equal( 43045383 %in% data_defaults@data.start, expected = FALSE )
})

# NOTE: This test indicates that destranding does not behave as we expect.
# If we want a minimum count of 10 after destranding, we need to set the
# minCount to 5 because updating after destranding occurs after filtering
# for minCount.
test_that('Test destranding rescue',{
  expect_equal( 43045074 %in% data_defaults@data.start, expected = TRUE )
  expect_equal( data_defaults@data.coverage[which(data_defaults@data.start == 43045074), 1], expected = 10)
})

test_that('Test destranding aggregation',{
  # Removal of sites moved to the other strand
  expect_equal( 43053298 %in% data_defaults@data.start, expected = FALSE )
  # Correct coverage sums
  expect_equal( data_defaults@data.coverage[which(data_defaults@data.start == 43053297), 1], expected = 96 )
  expect_equal( data_defaults@data.coverage[which(data_defaults@data.start == 43053297), 2], expected = 207 )
  # Correct numCs
  expect_equal( data_defaults@data.numCs[which(data_defaults@data.start == 43053297), 1], expected = 87 )
})

test_that('Test for correct nrows in methylSigData object',{
  expect_equal( nrow(data_defaults@data.coverage), expected = 8 )
})

################################################################################
# Test with SNP removal

data_snps = methylSigReadData(
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
  filterSNPs = TRUE,
  num.cores = 1,
  quiet = TRUE)

test_that('Test SNP removal',{
  expect_equal( 9419909 %in% data_snps@data.start, expected = FALSE )
})

test_that('Test for correct nrows in methylSigData object',{
  expect_equal( length(data_snps@data.start), expected = 7 )
})

################################################################################
# Test with strandedness (destranded = FALSE)
# And reduced minCount to 5

data_stranded = methylSigReadData(
  fileList = files,
  sample.ids = c('test_1','test_2'),
  assembly = 'hg19',
  pipeline = 'bismark and methylKit',
  header = FALSE,
  context = 'CpG',
  resolution = "base",
  treatment = c(1,0),
  destranded = FALSE,
  maxCount = 500,
  minCount = 5,
  filterSNPs = FALSE,
  num.cores = 1,
  quiet = TRUE)

test_that('Test for correct nrows in methylSigData object',{
  expect_equal( length(data_stranded@data.start), expected = 10)
})

################################################################################
# Test with destranded
# And maxCount lowered to 200

data_ceiling = methylSigReadData(
  fileList = files,
  sample.ids = c('test_1','test_2'),
  assembly = 'hg19',
  pipeline = 'bismark and methylKit',
  header = FALSE,
  context = 'CpG',
  resolution = "base",
  treatment = c(1,0),
  destranded = TRUE,
  maxCount = 200,
  minCount = 10,
  filterSNPs = FALSE,
  num.cores = 1,
  quiet = TRUE)

test_that('Test for removal of combined > maxCount',{
  expect_equal( 43053297 %in% data_ceiling@data.start, expected = FALSE)
})

test_that('Test for correct nrows in methylSigData object',{
  expect_equal( length(data_ceiling@data.start), expected = 7 )
})
