library(bsseq)
library(GenomicRanges)

dir.create('inst/extdata', recursive = TRUE)

########################################

#---------CG-------------CG-------------CG------------------C--------------CG
#---------GC-------------GC-------------GC------------------G--------------GC

gr1 = GRanges(
    seqnames = rep.int('chr1', 9),
    ranges = IRanges(
        start = c(10,11,25,26,40,41,60,75,76),
        end = c(11,12,26,27,41,42,61,76,77)),
    strand = c('+','-','+','-','+','-','+','+','-')
)

cov1 = matrix(
    data = c(5,5,30,70,10,20,40,1000,1500),
    ncol = 1
)

meth1 = matrix(
    data = c(4,4,0,5,9,19,35,900,1400),
    ncol = 1
)

########################################

#---------CG-------------CG-------------CG-------CG---------C--------------CG
#---------GC-------------GC-------------GC-------GC---------G--------------GC

gr2 = GRanges(
    seqnames = rep.int('chr1', 11),
    ranges = IRanges(
        start = c(10,11,25,26,40,41,50,51,60,75,76),
        end = c(11,12,26,27,41,42,51,52,61,76,77)),
    strand = c('+','-','+','-','+','-','+','-','+','+','-')
)

cov2 = matrix(
    data = c(10,10,50,50,15,35,5,5,20,100,200),
    ncol = 1
)

meth2 = matrix(
    data = c(9,9,1,5,14,34,5,5,15,99,199),
    ncol = 1
)

########################################

#---------C--------------C--------------C-------------------C--------------C-

gr3 = GRanges(
    seqnames = rep.int('chr1', 5),
    ranges = IRanges(
        start = c(10,25,40,60,75),
        end = c(12,27,42,62,77))
)

cov3 = matrix(
    data = c(10,100,30,40,2500),
    ncol = 1
)

meth3 = matrix(
    data = c(8,5,28,35,2300),
    ncol = 1
)

########################################

#---------C--------------C--------------C--------C----------C--------------C-

gr4 = GRanges(
    seqnames = rep.int('chr1', 6),
    ranges = IRanges(
        start = c(10,25,40,50,60,75),
        end = c(12,27,42,52,62,77))
)

cov4 = matrix(
    data = c(20,100,50,10,20,300),
    ncol = 1
)

meth4 = matrix(
    data = c(18,6,48,10,15,298),
    ncol = 1
)

########################################

bsseq1 = BSseq(
    M = meth1,
    Cov = cov1,
    gr = gr1,
    sampleNames = 'test1'
)

bsseq2 = BSseq(
    M = meth2,
    Cov = cov2,
    gr = gr2,
    sampleNames = 'test2'
)

bsseq3 = BSseq(
    M = meth3,
    Cov = cov3,
    gr = gr3,
    sampleNames = 'test3'
)

bsseq4 = BSseq(
    M = meth4,
    Cov = cov4,
    gr = gr4,
    sampleNames = 'test4'
)

#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 coverage
#         5              30             10       0          40             1000
#          5              70             20       0                         1500
# test2 coverage
#         10             50             15       5          20             100
#          10             50             35       5                         200
# test1 methylation
#         4              0              9        0          35             900
#          4              5              19       0                         1400
# test2 methylation
#         9              1              14       5          15             99
#          9              5              34       5                         199
bsseq_with_strand = combine(bsseq1, bsseq2)

#---------C--------------C--------------C--------C----------C--------------C-
# test1 coverage
#         10             100            30       0          40             2500
# test2 coverage
#         20             100            50       10         20             300
# test1 methylation
#         8              5              28       0          35             2300
# test2 methylation
#         18             6              48       10         15             298
bsseq_without_strand = combine(bsseq3, bsseq4)

########################################

df1 = data.frame(gr1)
df1$meth = as.numeric(meth1)
df1$unmeth = as.numeric(cov1 - meth1)
df1$perc = as.numeric(meth1 / cov1) * 100
df1$dicontext = rep.int('CG', 9)

df2 = data.frame(gr2)
df2$meth = as.numeric(meth2)
df2$unmeth = as.numeric(cov2 - meth2)
df2$perc = as.numeric(meth2 / cov2) * 100
df2$dicontext = rep.int('CG', 11)

df3 = data.frame(gr3)
df3$meth = as.numeric(meth3)
df3$unmeth = as.numeric(cov3 - meth3)
df3$perc = as.numeric(meth3 / cov3) * 100
df3$dicontext = rep.int('CG', 5)

df4 = data.frame(gr4)
df4$meth = as.numeric(meth4)
df4$unmeth = as.numeric(cov4 - meth4)
df4$perc = as.numeric(meth4 / cov4) * 100
df4$dicontext = rep.int('CG', 6)

########################################

cov_cols = c('seqnames','start','end','perc','meth','unmeth')

bis_cov_file1 = './inst/extdata/bis_cov1.cov'
bis_cov_file2 = './inst/extdata/bis_cov2.cov'
bis_cov_file3 = './inst/extdata/bis_cov3.cov'
bis_cov_file4 = './inst/extdata/bis_cov4.cov'

bis_cov1 = df1[, cov_cols]
bis_cov2 = df2[, cov_cols]
bis_cov3 = df3[, cov_cols]
bis_cov4 = df4[, cov_cols]

write.table(
    x = bis_cov1, file = bis_cov_file1,
    quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

write.table(
    x = bis_cov2, file = bis_cov_file2,
    quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

write.table(
    x = bis_cov3, file = bis_cov_file3,
    quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

write.table(
    x = bis_cov4, file = bis_cov_file4,
    quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

########################################

# report_cols = c('seqnames','start','strand','perc','meth','unmeth','dicontext')
#
# bis_report_file1 = './inst/extdata/bis_report1.txt'
# bis_report_file2 = './inst/extdata/bis_report2.txt'
#
# bis_report1 = df1[, report_cols]
# bis_report2 = df2[, report_cols]
#
# write.table(
#     x = bis_report1, file = bis_report_file1,
#     quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
#
# write.table(
#     x = bis_report2, file = bis_report_file2,
#     quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

########################################

bsseq_cov_with_strand = read.bismark(
    files = c(bis_cov_file1, bis_cov_file2),
    colData = data.frame(row.names = c('test1','test2')),
    rmZeroCov = FALSE,
    strandCollapse = FALSE
)

bsseq_cov_without_strand = read.bismark(
    files = c(bis_cov_file3, bis_cov_file4),
    colData = data.frame(row.names = c('test3','test4')),
    rmZeroCov = FALSE,
    strandCollapse = FALSE
)

# bsseq_report_with_strand = read.bismark(
#     files = c(bis_report_file1, bis_report_file2),
#     colData = data.frame(row.names = c('test1','test2')),
#     rmZeroCov = FALSE,
#     strandCollapse = FALSE
# )
#
# bsseq_report_without_strand = read.bismark(
#     files = c(bis_report_file1, bis_report_file2),
#     colData = data.frame(row.names = c('test1','test2')),
#     rmZeroCov = FALSE,
#     strandCollapse = TRUE
# )

usethis::use_data(
    bsseq_cov_with_strand,
    bsseq_cov_without_strand)

# usethis::use_data(
#     bsseq_report_with_strand,
#     bsseq_report_without_strand)

########################################

#----[---------]---------[--------------]----[---------]--------------[---------]
tile_regions_gr1 = GRanges(
    seqnames = c('chr1','chr1','chr1','chr1'),
    ranges = IRanges(
        start = c(5,25,45,70),
        end = c(15,40,55,80)
    )
)

### bsseq_cov_with_strand
## coverage
#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 coverage
#         5              30             10       0          40             1000
#          5              70             20       0                         1500
# test2 coverage
#         10             50             15       5          20             100
#          10             50             35       5                         200
#----[---------]---------[--------------]----[---------]--------------[---------]
# tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr1) coverage
#        10                     110             0                         2500
#        20                     115             10                        300

## methylation
#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 methylation
#         4              0              9        0          35             900
#          4              5              19       0                         1400
# test2 methylation
#         9              1              14       5          15             99
#          9              5              34       5                         199
#----[---------]---------[--------------]----[---------]--------------[---------]
# tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr1) methylation
#        8                     14               0                         2300
#        18                    20               10                        298

bsseq_with_strand_tile_gr1 = BSseq(
    gr = tile_regions_gr1,
    Cov = matrix(c(10,110,0,2500,20,115,10,300), ncol = 2),
    M = matrix(c(8,14,0,2300,18,20,10,298), ncol = 2),
    pData = data.frame(row.names = c('test1','test2')),
    sampleNames = c('test1','test2')
)

### bsseq_cov_without_strand
## coverage
#---------C--------------C--------------C--------C----------C--------------C-
# test1 coverage
#         10             100            30       0          40             2500
# test2 coverage
#         20             100            50       10         20             300
#----[---------]---------[--------------]----[---------]--------------[---------]
# tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr1)
#         10                    130               0                        2500
#         20                    150               10                       300

## methylation
#---------C--------------C--------------C--------C----------C--------------C-
# test1 methylation
#         8              5              28       0          35             2300
# test2 methylation
#         18             6              48       10         15             298
#----[---------]---------[--------------]----[---------]--------------[---------]
# tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr1)
#         8                    33               0                        2300
#         18                   54               10                       298

bsseq_without_strand_tile_gr1 = BSseq(
    gr = tile_regions_gr1,
    Cov = matrix(c(10,130,0,2500,20,150,10,300), ncol = 2),
    M = matrix(c(8,33,0,2300,18,54,10,298), ncol = 2),
    pData = data.frame(row.names = c('test3','test4')),
    sampleNames = c('test3','test4')
)

########################################
#----[------------------------]----[----------------------------------]----[----]
tile_regions_gr2 = GRanges(
    seqnames = c('chr1','chr1','chr1'),
    ranges = IRanges(
        start = c(5,35,75),
        end = c(30,70,80)
    )
)

### bsseq_cov_with_strand
## coverage
#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 coverage
#         5              30             10       0          40             1000
#          5              70             20       0                         1500
# test2 coverage
#         10             50             15       5          20             100
#          10             50             35       5                         200
#----[------------------------]----[----------------------------------]----[----]
# tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr2) coverage
#           110                                 70                          2500
#           120                                 80                          300

## methylation
#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 methylation
#         4              0              9        0          35             900
#          4              5              19       0                         1400
# test2 methylation
#         9              1              14       5          15             99
#          9              5              34       5                         199
#----[------------------------]----[----------------------------------]----[----]
# tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr2) methylation
#           13                                 63                            2300
#           24                                 73                            298
bsseq_with_strand_tile_gr2 = BSseq(
    gr = tile_regions_gr2,
    Cov = matrix(c(110,70,2500,120,80,300), ncol = 2),
    M = matrix(c(13,63,2300,24,73,298), ncol = 2),
    pData = data.frame(row.names = c('test1','test2')),
    sampleNames = c('test1','test2')
)

### bsseq_cov_without_strand
## coverage
#---------C--------------C--------------C--------C----------C--------------C-
# test1 coverage
#         10             100            30       0          40             2500
# test2 coverage
#         20             100            50       10         20             300
#----[------------------------]----[----------------------------------]----[----]
# tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr2)
#           110                                 70                          2500
#           120                                 80                          300

## methylation
#---------C--------------C--------------C--------C----------C--------------C-
# test1 methylation
#         8              5              28       0          35             2300
# test2 methylation
#         18             6              48       10         15             298
#----[------------------------]----[----------------------------------]----[----]
# tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr2)
#           13                                 63                            2300
#           24                                 73                            298
bsseq_without_strand_tile_gr2 = BSseq(
    gr = tile_regions_gr2,
    Cov = matrix(c(110,70,2500,120,80,300), ncol = 2),
    M = matrix(c(13,63,2300,24,73,298), ncol = 2),
    pData = data.frame(row.names = c('test3','test4')),
    sampleNames = c('test3','test4')
)

########################################
#----[--------------------------------------------------------------------------]
tile_regions_gr3 = GRanges(
    seqnames = c('chr1'),
    ranges = IRanges(
        start = c(5),
        end = c(80)
    )
)

### bsseq_cov_with_strand
## coverage
#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 coverage
#         5              30             10       0          40             1000
#          5              70             20       0                         1500
# test2 coverage
#         10             50             15       5          20             100
#          10             50             35       5                         200
#----[--------------------------------------------------------------------------]
# tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr3) coverage
#                                       2680
#                                       500

## methylation
#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 methylation
#         4              0              9        0          35             900
#          4              5              19       0                         1400
# test2 methylation
#         9              1              14       5          15             99
#          9              5              34       5                         199
#----[--------------------------------------------------------------------------]
# tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr3) methylation
#                                       2376
#                                       395
bsseq_with_strand_tile_gr3 = BSseq(
    gr = tile_regions_gr3,
    Cov = matrix(c(2680,500), ncol = 2),
    M = matrix(c(2376,395), ncol = 2),
    pData = data.frame(row.names = c('test1','test2')),
    sampleNames = c('test1','test2')
)

### bsseq_cov_without_strand
## coverage
#---------C--------------C--------------C--------C----------C--------------C-
# test1 coverage
#         10             100            30       0          40             2500
# test2 coverage
#         20             100            50       10         20             300
#----[--------------------------------------------------------------------------]
# tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr3)
#                                       2680
#                                       500

## methylation
#---------C--------------C--------------C--------C----------C--------------C-
# test1 methylation
#         8              5              28       0          35             2300
# test2 methylation
#         18             6              48       10         15             298
#----[--------------------------------------------------------------------------]
# tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr3)
#                                       2376
#                                       395
bsseq_without_strand_tile_gr3 = BSseq(
    gr = tile_regions_gr3,
    Cov = matrix(c(2680,500), ncol = 2),
    M = matrix(c(2376,395), ncol = 2),
    pData = data.frame(row.names = c('test3','test4')),
    sampleNames = c('test3','test4')
)

########################################
#----[--------------]------------------------------------------------------------
tile_regions_gr4 = GRanges(
    seqnames = c('chr1'),
    ranges = IRanges(
        start = c(5),
        end = c(20)
    )
)

### bsseq_cov_with_strand
## coverage
#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 coverage
#         5              30             10       0          40             1000
#          5              70             20       0                         1500
# test2 coverage
#         10             50             15       5          20             100
#          10             50             35       5                         200
#----[--------------]------------------------------------------------------------
# tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr4) coverage
#          10
#          20

## methylation
#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 methylation
#         4              0              9        0          35             900
#          4              5              19       0                         1400
# test2 methylation
#         9              1              14       5          15             99
#          9              5              34       5                         199
#----[--------------]------------------------------------------------------------
# tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr4) methylation
#          8
#          18
bsseq_with_strand_tile_gr4 = BSseq(
    gr = tile_regions_gr4,
    Cov = matrix(c(10,20), ncol = 2),
    M = matrix(c(8,18), ncol = 2),
    pData = data.frame(row.names = c('test1','test2')),
    sampleNames = c('test1','test2')
)

### bsseq_cov_without_strand
## coverage
#---------C--------------C--------------C--------C----------C--------------C-
# test1 coverage
#         10             100            30       0          40             2500
# test2 coverage
#         20             100            50       10         20             300
#----[--------------]------------------------------------------------------------
# tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr4)
#          10
#          20

## methylation
#---------C--------------C--------------C--------C----------C--------------C-
# test1 methylation
#         8              5              28       0          35             2300
# test2 methylation
#         18             6              48       10         15             298
#----[--------------]------------------------------------------------------------
# tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr4)
#          8
#          18
bsseq_without_strand_tile_gr4 = BSseq(
    gr = tile_regions_gr4,
    Cov = matrix(c(10,20), ncol = 2),
    M = matrix(c(8,18), ncol = 2),
    pData = data.frame(row.names = c('test3','test4')),
    sampleNames = c('test3','test4')
)

########################################
#----[---]-----------------------------------------------------------------------
tile_regions_gr5 = GRanges(
    seqnames = c('chr1'),
    ranges = IRanges(
        start = c(5),
        end = c(9)
    )
)

### bsseq_cov_with_strand
## coverage
#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 coverage
#         5              30             10       0          40             1000
#          5              70             20       0                         1500
# test2 coverage
#         10             50             15       5          20             100
#          10             50             35       5                         200
#----[---]-----------------------------------------------------------------------
# tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr5) coverage
#
#

## methylation
#---------CG-------------CG-------------CG-------CG---------C--------------CG
# test1 methylation
#         4              0              9        0          35             900
#          4              5              19       0                         1400
# test2 methylation
#         9              1              14       5          15             99
#          9              5              34       5                         199
#----[---]-----------------------------------------------------------------------
# tile_by_regions(bs = bsseq_cov_with_strand, gr = tile_regions_gr5) methylation
#
#
bsseq_with_strand_tile_gr5 = BSseq(
    gr = tile_regions_gr5,
    Cov = matrix(c(0,0), ncol = 2),
    M = matrix(c(0,0), ncol = 2),
    pData = data.frame(row.names = c('test1','test2')),
    sampleNames = c('test1','test2')
)

### bsseq_cov_without_strand
## coverage
#---------C--------------C--------------C--------C----------C--------------C-
# test1 coverage
#         10             100            30       0          40             2500
# test2 coverage
#         20             100            50       10         20             300
#----[---]-----------------------------------------------------------------------
# tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr5)
#
#

## methylation
#---------C--------------C--------------C--------C----------C--------------C-
# test1 methylation
#         8              5              28       0          35             2300
# test2 methylation
#         18             6              48       10         15             298
#----[---]-----------------------------------------------------------------------
# tile_by_regions(bs = bsseq_cov_without_strand, gr = tile_regions_gr5)
#
#
bsseq_without_strand_tile_gr5 = BSseq(
    gr = tile_regions_gr5,
    Cov = matrix(c(0,0), ncol = 2),
    M = matrix(c(0,0), ncol = 2),
    pData = data.frame(row.names = c('test3','test4')),
    sampleNames = c('test3','test4')
)

########################################

usethis::use_data(
    tile_regions_gr1,
    tile_regions_gr2,
    tile_regions_gr3,
    tile_regions_gr4,
    tile_regions_gr5,
    bsseq_with_strand_tile_gr1,
    bsseq_without_strand_tile_gr1,
    bsseq_with_strand_tile_gr2,
    bsseq_without_strand_tile_gr2,
    bsseq_with_strand_tile_gr3,
    bsseq_without_strand_tile_gr3,
    bsseq_with_strand_tile_gr4,
    bsseq_without_strand_tile_gr4,
    bsseq_with_strand_tile_gr5,
    bsseq_without_strand_tile_gr5,
    internal = TRUE,
    overwrite = TRUE)
