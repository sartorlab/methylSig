library(bsseq)

bis_cov_file1 = './inst/extdata/bis_cov1.cov'
bis_cov_file2 = './inst/extdata/bis_cov2.cov'
bis_cov_file3 = './inst/extdata/bis_cov3.cov'
bis_cov_file4 = './inst/extdata/bis_cov4.cov'

#---------CG-------------CG-------------CG--------CG--------C--------------CG
# test1 coverage
#         5              30             10        0         40             1000
#          5              70             20        0                         1500
# test2 coverage
#         10             50             15        5         20             100
#          10             50             35        5                         200
# test1 methylation
#         4              0              9         0         35             900
#          4              5              19        0                         1400
# test2 methylation
#         9              1              14        5         15             99
#          9              5              34        5                         199
bsseq_stranded = read.bismark(
    files = c(bis_cov_file1, bis_cov_file2),
    colData = data.frame(row.names = c('test1','test2')),
    rmZeroCov = FALSE,
    strandCollapse = FALSE
)

#---------C--------------C--------------C---------C---------C--------------C-
# test1 coverage
#         10             100            30        0         40             2500
# test2 coverage
#         20             100            50        10        20             300
# test1 methylation
#         8              5              28        0         35             2300
# test2 methylation
#         18             6              48        10        15             298
bsseq_destranded = read.bismark(
    files = c(bis_cov_file3, bis_cov_file4),
    colData = data.frame(row.names = c('test3','test4')),
    rmZeroCov = FALSE,
    strandCollapse = FALSE
)

usethis::use_data(
    bsseq_stranded,
    bsseq_destranded)
