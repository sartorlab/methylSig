# The index of C > T SNPs will be a GenomicRanges object for efficient overlap finding
library(GenomicRanges)

# Annotate with hg19 chromosome lengths
chr_lengths = read.table('~/latte/Methylation/Data/chromInfo_hg19.txt', header=F, sep='\t', stringsAsFactors=F)

# Grab the vcf, but only consider the first two columns
setwd('~/latte/vcf')
vcf = read.table('filtered_AF_0.05_CT_SNPs.vcf',header=F,sep='\t',quote='',comment.char='', colClasses=c('character','numeric','NULL','NULL','NULL','NULL','NULL','NULL'),stringsAsFactors=F)
colnames(vcf) = c('chrom','start')
vcf$chrom = paste('chr',vcf$chrom,sep='')

CT_SNPs_hg19 = GRanges(seqnames=vcf$chrom, ranges=IRanges(start=vcf$start, end=vcf$start))
seqlengths(CT_SNPs_hg19) = chr_lengths[match(names(seqlengths(CT_SNPs_hg19)), chr_lengths[,1]), 2]
genome(CT_SNPs_hg19) = 'hg19'

save(CT_SNPs_hg19, file='CT_SNPs_hg19.RData')
write.table(vcf, file='CT_SNPs_hg19.txt', sep='\t', row.names=F, quote=F)
