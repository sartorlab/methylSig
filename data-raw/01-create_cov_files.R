library(bsseq)
library(GenomicRanges)

dir.create('inst/extdata', recursive = TRUE)

################################################################################

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

################################################################################

#---------CG-------------CG-------------CG--------CG--------C--------------CG
#---------GC-------------GC-------------GC--------GC--------G--------------GC

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

################################################################################

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

################################################################################

#---------C--------------C--------------C---------C---------C--------------C-

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

################################################################################

df1 = data.frame(gr1)
df1$meth = as.numeric(meth1)
df1$unmeth = as.numeric(cov1 - meth1)
df1$perc = as.numeric(meth1 / cov1) * 100

df2 = data.frame(gr2)
df2$meth = as.numeric(meth2)
df2$unmeth = as.numeric(cov2 - meth2)
df2$perc = as.numeric(meth2 / cov2) * 100

df3 = data.frame(gr3)
df3$meth = as.numeric(meth3)
df3$unmeth = as.numeric(cov3 - meth3)
df3$perc = as.numeric(meth3 / cov3) * 100

df4 = data.frame(gr4)
df4$meth = as.numeric(meth4)
df4$unmeth = as.numeric(cov4 - meth4)
df4$perc = as.numeric(meth4 / cov4) * 100

################################################################################

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
