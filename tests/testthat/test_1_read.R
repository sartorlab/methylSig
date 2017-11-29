context('Test methylSigReadData')

# Lay out what you expect to happen in all cases

########################
### test_1.txt
# destrand = T, filterSNPs = T, min/max = 10/500: cov = -; M = - # at 9413839
# destrand = T, filterSNPs = F, min/max = 10/500: cov = 45; M = 40 # at 9413839
# destrand = F, filterSNPs = T, min/max = 10/500: cov = -/25; M = -/20
# destrand = F, filterSNPs = F, min/max = 10/500: cov = 20/25; M = 0/5
# 9413839 is a CT SNP
# chr21	9413839	+	20	0	CG	CGG
# chr21	9413840	-	20	5	CG	CGG

# destrand = T, filterSNPs = T, min/max = 10/500: cov = 34; M = 34 # at 9419909
# destrand = T, filterSNPs = F, min/max = 10/500: cov = 34; M = 34 # at 9419909
# destrand = F, filterSNPs = T, min/max = 10/500: cov = 14/-; M = 14/-
# destrand = F, filterSNPs = F, min/max = 10/500: cov = 14/20; M = 14/20
# 9419909 is a CT SNP
# chr21	9419908	+	14	0	CG	CGG
# chr21	9419909	-	20	0	CG	CGG

# destrand = T, filterSNPs = T, min/max = 10/500: cov = 0; M = 0 (kept because test_2.txt has identical location that is kept)
# destrand = T, filterSNPs = F, min/max = 10/500: cov = 0; M = 0 (kept because test_2.txt has identical location that is kept)
# destrand = F, filterSNPs = T, min/max = 10/500: cov = 0; M = 0 (kept because test_2.txt has identical location that is kept)
# destrand = F, filterSNPs = F, min/max = 10/500: cov = 0; M = 0 (kept because test_2.txt has identical location that is kept)
# No SNP
# chr21	43052356	-	550	0	CG	CGG

# destrand = T, filterSNPs = T, min/max = 10/500: cov = 96; M = 87 # at 43053297
# destrand = T, filterSNPs = F, min/max = 10/500: cov = 96; M = 87 # at 43053297
# destrand = F, filterSNPs = T, min/max = 10/500: cov = 24/72; M = 23/64
# destrand = F, filterSNPs = F, min/max = 10/500: cov = 24/72; M = 23/64
# No SNP
# chr21	43053297	+	23	1	CG	CGG
# chr21	43053298	-	64	8	CG	CGG

# destrand = T, filterSNPs = T, min/max = 10/500: cov = -; M = - (removed because test_2.txt has identical and 0 in both)
# destrand = T, filterSNPs = F, min/max = 10/500: cov = -; M = - (removed because test_2.txt has identical and 0 in both)
# destrand = F, filterSNPs = T, min/max = 10/500: cov = -; M = - (removed because test_2.txt has identical and 0 in both)
# destrand = F, filterSNPs = F, min/max = 10/500: cov = -; M = - (removed because test_2.txt has identical and 0 in both)
# No SNP
# chr21	43053323	+	4	0	CG	CGG

########################
### test_2.txt
# destrand = T, filterSNPs = T, min/max = 10/500: cov = -; M = - # at 9413839
# destrand = T, filterSNPs = F, min/max = 10/500: cov = 0; M = 0 # at 9413839
# destrand = F, filterSNPs = T, min/max = 10/500: cov = -/0; M = -/0
# destrand = F, filterSNPs = F, min/max = 10/500: cov = 0/0; M = 0/0
# EXCLUDE THESE FROM THIS FILE
# chr21	9413839	+	20	0	CG	CGG
# chr21	9413840	-	20	5	CG	CGG

# destrand = T, filterSNPs = T, min/max = 10/500: cov = 10; M = 10 # at 9419909
# destrand = T, filterSNPs = F, min/max = 10/500: cov = 10; M = 10  # at 9419909
# destrand = F, filterSNPs = T, min/max = 10/500: cov = 0/-; M = 0/-
# destrand = F, filterSNPs = F, min/max = 10/500: cov = 0/0; M = 0/0
# 9419910 is a CT SNP
# chr21	9419908	+	5	0	CG	CGG
# chr21	9419909	-	5	0	CG	CGG

# destrand = T, filterSNPs = T, min/max = 10/500: cov = 300; M = 300
# destrand = T, filterSNPs = F, min/max = 10/500: cov = 300; M = 300
# destrand = F, filterSNPs = T, min/max = 10/500: cov = 300; M = 300
# destrand = F, filterSNPs = F, min/max = 10/500: cov = 300; M = 300
# No SNP
# chr21	43052356	-	300	0	CG	CGG

# destrand = T, filterSNPs = T, min/max = 10/500: cov = 458; M = 434 # at 43053297 (kept because test_1.txt has identical location that is kept)
# destrand = T, filterSNPs = F, min/max = 10/500: cov = 458; M = 434 # at 43053297 (kept because test_1.txt has identical location that is kept)
# destrand = F, filterSNPs = T, min/max = 10/500: cov = 318/140; M = 300/134
# destrand = F, filterSNPs = F, min/max = 10/500: cov = 318/140; M = 300/134
# No SNP
# chr21	43053297	+	300	18	CG	CGG
# chr21	43053298	-	134	6	CG	CGG

# destrand = T, filterSNPs = T, min/max = 10/500: cov = -; M = - (removed because test_2.txt has identical and 0 in both)
# destrand = T, filterSNPs = F, min/max = 10/500: cov = -; M = - (removed because test_2.txt has identical and 0 in both)
# destrand = F, filterSNPs = T, min/max = 10/500: cov = -; M = - (removed because test_2.txt has identical and 0 in both)
# destrand = F, filterSNPs = F, min/max = 10/500: cov = -; M = - (removed because test_2.txt has identical and 0 in both)
# No SNP
# chr21	43053323	+	4	0	CG	CGG

result1_cov = matrix(c(34,0,96,10,300,458), nrow = 3, ncol = 2, byrow = FALSE)
result1_M = matrix(c(34,0,87,10,300,434), nrow = 3, ncol = 2, byrow = FALSE)

result2_cov = matrix(c(45,34,0,96,0,10,300,458), nrow = 4, ncol = 2, byrow = FALSE)
result2_M = matrix(c(40,34,0,87,0,10,300,434), nrow = 4, ncol = 2, byrow = FALSE)

result3_cov = matrix(c(25,14,0,24,72,0,0,300,318,140), nrow = 5, ncol = 2, byrow = FALSE)
result3_M = matrix(c(20,14,0,23,64,0,0,300,300,134), nrow = 5, ncol = 2, byrow = FALSE)

result4_cov = matrix(c(20,25,14,20,0,24,72,0,0,0,0,300,318,140), nrow = 7, ncol = 2, byrow = FALSE)
result4_M = matrix(c(20,20,14,20,0,23,64,0,0,0,0,300,300,134), nrow = 7, ncol = 2, byrow = FALSE)


################################################################################
# Test min/max filters

files = c(system.file('extdata', 'test_1.txt', package='methylSig'),
    system.file('extdata', 'test_2.txt', package='methylSig'))

test_that('data1 coverage and M matrices are as expected', {
    data1 = methylSigReadData(
        fileList = files,
        sample.ids = c('test_1', 'test_2'),
        assembly = 'hg19',
        destranded = TRUE,
        maxCount = 500,
        minCount = 10,
        filterSNPs = TRUE,
        num.cores = 1,
        fileType = 'cytosineReport')

    expect_true(all(as.matrix(bsseq::getCoverage(data1, type = 'Cov')) == result1_cov))
    expect_true(all(as.matrix(bsseq::getCoverage(data1, type = 'M')) == result1_M))
})

test_that('data2 coverage and M matrices are as expected', {
    data2 = methylSigReadData(
        fileList = files,
        sample.ids = c('test_1', 'test_2'),
        assembly = 'hg19',
        destranded = TRUE,
        maxCount = 500,
        minCount = 10,
        filterSNPs = FALSE,
        num.cores = 1,
        fileType = 'cytosineReport')

    expect_true(all(as.matrix(bsseq::getCoverage(data2, type = 'Cov')) == result2_cov))
    expect_true(all(as.matrix(bsseq::getCoverage(data2, type = 'M')) == result2_M))
})

test_that('data3 coverage and M matrices are as expected', {
    data3 = methylSigReadData(
        fileList = files,
        sample.ids = c('test_1', 'test_2'),
        assembly = 'hg19',
        destranded = FALSE,
        maxCount = 500,
        minCount = 10,
        filterSNPs = TRUE,
        num.cores = 1,
        fileType = 'cytosineReport')

    expect_true(all(as.matrix(bsseq::getCoverage(data3, type = 'Cov')) == result3_cov))
    expect_true(all(as.matrix(bsseq::getCoverage(data3, type = 'M')) == result3_M))
})

test_that('data3 coverage and M matrices are as expected', {
    data4 = methylSigReadData(
        fileList = files,
        sample.ids = c('test_1', 'test_2'),
        assembly = 'hg19',
        destranded = FALSE,
        maxCount = 500,
        minCount = 10,
        filterSNPs = FALSE,
        num.cores = 1,
        fileType = 'cytosineReport')

    expect_true(all(as.matrix(bsseq::getCoverage(data4, type = 'Cov')) == result4_cov))
    expect_true(all(as.matrix(bsseq::getCoverage(data4, type = 'M')) == result4_M))
})
