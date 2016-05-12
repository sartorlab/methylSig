context('Test methylSigReadData')

################################################################################
# Test min/max filters

files = system.file('extdata', 'test_minmax_1.txt', package='methylSig')

data1 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = FALSE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test minCount removal (data1)', {
	expect_equal( length(data1@data.start), expected = 5 )
	expect_equal( 3 %in% data1@data.start, expected = FALSE )
})

data2 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = TRUE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test minCount rescue (data2)', {
	expect_equal( length(data2@data.start), expected = 3 )
	expect_equal( 3 %in% data2@data.start, expected = TRUE )
	expect_equal( all( c(3,6,8) %in% data2@data.start ), expected = TRUE )
})

##############

files = system.file('extdata', 'test_minmax_2.txt', package='methylSig')

data3 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = FALSE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test maxCount removal (data3)', {
	expect_equal( length(data3@data.start), expected = 5 )
	expect_equal( 8 %in% data3@data.start, expected = FALSE )
})

data4 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = TRUE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test maxCount removal (data4)', {
	expect_equal( length(data4@data.start), expected = 2 )
	expect_equal( all( c(3,6) %in% data4@data.start ), expected = TRUE )
	expect_equal( 8 %in% data4@data.start, expected = FALSE )
	expect_equal( all( c(4,7,9) %in% data4@data.start ), expected = FALSE )
})

##############

files = c(
	system.file('extdata', 'test_minmax_1.txt', package='methylSig'),
	system.file('extdata', 'test_minmax_2.txt', package='methylSig')
	)

data5 = methylSigReadData(
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
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test minCount/maxCount site preservation (data5)', {
	expect_equal( length(data5@data.start), expected = 6 )
	expect_equal( data5@data.coverage[which(data5@data.start == 3), 1], expected = 0 )
	expect_equal( data5@data.coverage[which(data5@data.start == 3), 2], expected = 10 )
	expect_equal( data5@data.coverage[which(data5@data.start == 8), 1], expected = 200 )
	expect_equal( data5@data.coverage[which(data5@data.start == 8), 2], expected = 0 )
})

data6 = methylSigReadData(
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

test_that('Test minCount/maxCount site preservation (data6)', {
	expect_equal( length(data6@data.start), expected = 3 )
	expect_equal( data6@data.coverage[which(data6@data.start == 3), 1], expected = 19 )
	expect_equal( data6@data.coverage[which(data6@data.start == 3), 2], expected = 24 )
	expect_equal( data6@data.coverage[which(data6@data.start == 8), 1], expected = 400 )
	expect_equal( data6@data.coverage[which(data6@data.start == 8), 2], expected = 0 )
})

################################################################################
# Test %C + %T < 95%

files = system.file('extdata', 'test_freq_1.txt', package='methylSig')

data7 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = FALSE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test %C + %T < 95% removal (data7)', {
	expect_equal( length(data7@data.start), expected = 5 )
	expect_equal( 3 %in% data7@data.start, expected = FALSE )
})

data8 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = TRUE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

# Is this the behavior we want? In this case, destranding causes the data
# on the forward strand to be destroyed and shifts the data on the reverse
# strand to that position. Or do we want to do this check after destranding?
# That would cause the CpG at 3 to be lost completely.
test_that('Test %C + %T < 95% removal (data8)', {
	expect_equal( length(data8@data.start), expected = 3 )
	expect_equal( 3 %in% data8@data.start, expected = TRUE )
	expect_equal( data8@data.coverage[which(data8@data.start == 3), 1], expected = 14 )
})

##############

files = system.file('extdata', 'test_freq_2.txt', package='methylSig')

data9 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = FALSE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test %C + %T < 95% removal (data9)', {
	expect_equal( length(data9@data.start), expected = 5 )
	expect_equal( 9 %in% data9@data.start, expected = FALSE )
})

data10 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = TRUE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test %C + %T < 95% removal (data10)', {
	expect_equal( length(data10@data.start), expected = 3 )
	expect_equal( 8 %in% data10@data.start, expected = TRUE )
	expect_equal( data10@data.coverage[which(data10@data.start == 8), 1], expected = 200 )
})

##############

files = c(
	system.file('extdata', 'test_freq_1.txt', package='methylSig'),
	system.file('extdata', 'test_freq_2.txt', package='methylSig')
	)

data11 = methylSigReadData(
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
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test %C + %T < 95% removal (data11)', {
	expect_equal( length(data11@data.start), expected = 6 )
	expect_equal( data11@data.coverage[which(data11@data.start == 3), 1], expected = 0 )
	expect_equal( data11@data.coverage[which(data11@data.start == 3), 2], expected = 10 )
	expect_equal( data11@data.coverage[which(data11@data.start == 9), 1], expected = 200 )
	expect_equal( data11@data.coverage[which(data11@data.start == 9), 2], expected = 0 )
})

data12 = methylSigReadData(
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

test_that('Test %C + %T < 95% removal (data12)', {
	expect_equal( length(data12@data.start), expected = 3 )
	expect_equal( data12@data.coverage[which(data12@data.start == 3), 1], expected = 14 )
	expect_equal( data12@data.coverage[which(data12@data.start == 3), 2], expected = 24 )
	expect_equal( data12@data.coverage[which(data12@data.start == 8), 1], expected = 400 )
	expect_equal( data12@data.coverage[which(data12@data.start == 8), 2], expected = 200 )
})

################################################################################
# Test combinations

files = system.file('extdata', 'test_combo_1.txt', package='methylSig')

data13 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = FALSE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test combination (data13)', {
	expect_equal( length(data13@data.start), expected = 7 )
	expect_equal( all( c(1,3,4,8,10,30,31) %in% data13@data.start ), expected = FALSE )
	expect_equal( all( c(6,7,15,16,20,25,9411410) %in% data13@data.start ), expected = TRUE )
})

data14 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = FALSE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = TRUE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test combination (data14)', {
	expect_equal( length(data14@data.start), expected = 6 )
	expect_equal( all( c(1,3,4,8,10,30,31) %in% data14@data.start ), expected = FALSE )
	expect_equal( all( c(6,7,15,16,20,25) %in% data14@data.start ), expected = TRUE )
})

data15 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = TRUE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test combination (data15)', {
	expect_equal( length(data15@data.start), expected = 5 )
	expect_equal( all( c(1,4,8,10,15,25,30,31) %in% data15@data.start ), expected = FALSE )
	expect_equal( all( c(3,6,20,24,9411410) %in% data15@data.start ), expected = TRUE )
})

data16 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = TRUE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = TRUE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test combination (data16)', {
	expect_equal( length(data16@data.start), expected = 4 )
	expect_equal( all( c(1,4,8,10,15,25,30,31) %in% data16@data.start ), expected = FALSE )
	expect_equal( all( c(3,6,20,24) %in% data16@data.start ), expected = TRUE )
})

##############

files = system.file('extdata', 'test_combo_2.txt', package='methylSig')

data17 = methylSigReadData(
	fileList = files,
	sample.ids = c('test_1'),
	assembly = 'hg19',
	pipeline = 'bismark and methylKit',
	header = FALSE,
	context = 'CpG',
	resolution = "base",
	treatment = c(1),
	destranded = TRUE,
	maxCount = 500,
	minCount = 10,
	filterSNPs = TRUE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test combination (data17)', {
	expect_equal( length(data17@data.start), expected = 6 )
	expect_equal( all( c(10,16,25,30,31,48086569,48086570) %in% data17@data.start ), expected = FALSE )
	expect_equal( all( c(1,6,8,15,20,24) %in% data17@data.start ), expected = TRUE )
	expect_equal( data17@data.coverage[which(data17@data.start == 6), 1], expected = 55 )
	expect_equal( data17@data.coverage[which(data17@data.start == 15), 1], expected = 475 )
})

##############

files = c(
	system.file('extdata', 'test_combo_1.txt', package='methylSig'),
	system.file('extdata', 'test_combo_2.txt', package='methylSig')
	)

data18 = methylSigReadData(
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
	minCount = 10,
	filterSNPs = FALSE,
	num.cores = 1,
	quiet = TRUE)

test_that('Test combination (data18)', {
	expect_equal( length(data18@data.start), expected = 10 )
	expect_equal( all( c(10,30,31) %in% data18@data.start ), expected = FALSE )
	expect_equal( all( c(1,6,7,8,15,16,20,25,9411410,48086570) %in% data18@data.start ), expected = TRUE )
	expect_equal( data18@data.coverage[which(data18@data.start == 1), 1], expected = 0 )
	expect_equal( data18@data.coverage[which(data18@data.start == 1), 2], expected = 14 )
})

data19 = methylSigReadData(
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

test_that('Test combination (data19)', {
	expect_equal( length(data19@data.start), expected = 9 )
	expect_equal( all( c(4,7,10,16,25,30,31,48086570) %in% data19@data.start ), expected = FALSE )
	expect_equal( all( c(1,3,6,8,15,20,24,9411410,48086569) %in% data19@data.start ), expected = TRUE )
	expect_equal( data19@data.coverage[which(data19@data.start == 3), 1], expected = 12 )
	expect_equal( data19@data.coverage[which(data19@data.start == 3), 2], expected = 0 )
	expect_equal( data19@data.coverage[which(data19@data.start == 6), 1], expected = 36 )
	expect_equal( data19@data.coverage[which(data19@data.start == 6), 2], expected = 55 )
	expect_equal( data19@data.coverage[which(data19@data.start == 8), 1], expected = 0 )
	expect_equal( data19@data.coverage[which(data19@data.start == 8), 2], expected = 400 )
	expect_equal( data19@data.coverage[which(data19@data.start == 15), 1], expected = 0 )
	expect_equal( data19@data.coverage[which(data19@data.start == 15), 2], expected = 475 )
	expect_equal( data19@data.coverage[which(data19@data.start == 20), 1], expected = 30 )
	expect_equal( data19@data.coverage[which(data19@data.start == 20), 2], expected = 80 )
	expect_equal( data19@data.coverage[which(data19@data.start == 24), 1], expected = 50 )
	expect_equal( data19@data.coverage[which(data19@data.start == 24), 2], expected = 25 )
	expect_equal( data19@data.coverage[which(data19@data.start == 9411410), 1], expected = 72 )
	expect_equal( data19@data.coverage[which(data19@data.start == 9411410), 2], expected = 0 )
	expect_equal( data19@data.coverage[which(data19@data.start == 48086569), 1], expected = 0 )
	expect_equal( data19@data.coverage[which(data19@data.start == 48086569), 2], expected = 34 )
})

data20 = methylSigReadData(
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

test_that('Test combination (data20)', {
	expect_equal( length(data20@data.start), expected = 7 )
	expect_equal( all( c(4,7,10,16,25,30,31,9411410,48086569,48086570) %in% data20@data.start ), expected = FALSE )
	expect_equal( all( c(1,3,6,8,15,20,24) %in% data20@data.start ), expected = TRUE )
	expect_equal( data20@data.coverage[which(data20@data.start == 3), 1], expected = 12 )
	expect_equal( data20@data.coverage[which(data20@data.start == 3), 2], expected = 0 )
	expect_equal( data20@data.coverage[which(data20@data.start == 6), 1], expected = 36 )
	expect_equal( data20@data.coverage[which(data20@data.start == 6), 2], expected = 55 )
	expect_equal( data20@data.coverage[which(data20@data.start == 8), 1], expected = 0 )
	expect_equal( data20@data.coverage[which(data20@data.start == 8), 2], expected = 400 )
	expect_equal( data20@data.coverage[which(data20@data.start == 15), 1], expected = 0 )
	expect_equal( data20@data.coverage[which(data20@data.start == 15), 2], expected = 475 )
	expect_equal( data20@data.coverage[which(data20@data.start == 20), 1], expected = 30 )
	expect_equal( data20@data.coverage[which(data20@data.start == 20), 2], expected = 80 )
	expect_equal( data20@data.coverage[which(data20@data.start == 24), 1], expected = 50 )
	expect_equal( data20@data.coverage[which(data20@data.start == 24), 2], expected = 25 )
})
