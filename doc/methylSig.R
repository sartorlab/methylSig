## ---- echo=FALSE-----------------------------------------------------------
library(methylSig)
library(rtracklayer)

## ---- eval=FALSE-----------------------------------------------------------
#  devtools::install_github('sartorlab/methylSig')

## --------------------------------------------------------------------------
# The following bismark cytosine reports are included in inst/extdata
files = c(
    system.file('extdata', 'MDAMB_231_1DR.txt.gz', package='methylSig'),
    system.file('extdata', 'MDAMB_231_1DS.txt.gz', package='methylSig'),
    system.file('extdata', 'MDAMB_231_2DR.txt.gz', package='methylSig'),
    system.file('extdata', 'MDAMB_231_2DS.txt.gz', package='methylSig'),
    system.file('extdata', 'MDAMB_231_3DR.txt.gz', package='methylSig'),
    system.file('extdata', 'MDAMB_231_3DS.txt.gz', package='methylSig'))

sample.ids = basename(files)
sample.ids = gsub('.txt.gz', '', sample.ids)

# Build a pData matrix with columns for the samples, group memberships, and phenotype data
pData = data.frame(
    Sample_Names = sample.ids,
    DR_vs_DS = relevel(factor(c('DR','DS','DR','DS','DR','DS')), ref = 'DS'),
    row.names = sample.ids,
    stringsAsFactors = FALSE)

meth = methylSigReadData(
    fileList = files,
    pData = pData,
    assembly = 'hg19',
    destranded = TRUE,
    maxCount = 500,
    minCount = 10,
    filterSNPs = TRUE,
    num.cores = 1,
    fileType = 'cytosineReport')

print(meth)

## --------------------------------------------------------------------------
### Test on CpGs
result = methylSigCalc(
    meth = meth,
    comparison = 'DR_vs_DS',
    dispersion = 'both',
    local.info = FALSE,
    local.winsize = 200,
    min.per.group = c(3,3),
    weightFunc = methylSig_weightFunc,
    T.approx = TRUE,
    num.cores = 1)

print(result)

## --------------------------------------------------------------------------
# Must create a design matrix
design1 = data.frame(group = bsseq::pData(meth)$DR_vs_DS)

print(design1)

# NOTE this model has an intercept
contrast_intercept = matrix(c(0,1), ncol = 1)
result_dss_intercept = methylSigDSS(
    meth = meth,
    design = design1,
    formula = '~ group',
    contrast = contrast_intercept,
    group.term = 'group',
    min.per.group=c(3,3))

print(result_dss_intercept)

## --------------------------------------------------------------------------
# Add a covariate column, note specification as a factor, but can
# also use a numeric covariate
design2 = data.frame(
    group = bsseq::pData(meth)$DR_vs_DS,
    subject = factor(c(1,1,2,2,3,3)))

print(design2)

# NOTE the contrast vector has as many entries as the sum of the
# levels in group and subject, in the formula.
contrast_covariates = matrix(c(0,1,0,0), ncol = 1)
result_dss_covariates = methylSigDSS(
    meth = meth,
    design = design2,
    formula = '~ group + subject',
    contrast = contrast_covariates,
    group.term = 'group',
    min.per.group=c(3,3))

print(result_dss_covariates)

## --------------------------------------------------------------------------
### Test on 10000bp windows
windowed_meth = methylSigTile(meth, tiles = NULL, win.size = 10000)

tiled_result = methylSigCalc(
    meth = windowed_meth,
    comparison = 'DR_vs_DS',
    dispersion = 'both',
    local.info = FALSE,
    local.winsize = 200,
    min.per.group = c(3,3),
    weightFunc = methylSig_weightFunc,
    T.approx = TRUE,
    num.cores = 1)

print(tiled_result)

## --------------------------------------------------------------------------
### Test on CpG islands
library(annotatr)

cpg_islands = annotatr::build_annotations(genome = 'hg19', annotations = 'hg19_cpg_islands')

cpg_island_meth = methylSigTile(meth, tiles = cpg_islands)

cpg_island_result = methylSigCalc(
    meth = cpg_island_meth,
    comparison = 'DR_vs_DS',
    dispersion = 'both',
    local.info = FALSE,
    local.winsize = 200,
    min.per.group = c(3,3),
    weightFunc = methylSig_weightFunc,
    T.approx = TRUE,
    num.cores = 1)

print(cpg_island_result)

## --------------------------------------------------------------------------
# Get CpG island annotations from built-in data they could be built with the following:
# cpg_annots = annotatr::build_annotations(genome = 'hg19', annotations = c('hg19_cpg_islands'))
utils::data(sample_data, package = 'methylSig')

# Determine what CpGs should be considered significant
dmcList = result$fdr < 0.05 & abs(result$meth.diff) > 25

annotated_result = methylSigAnnotation(myDiff = result, dmcList = dmcList, annotations = cpg_annots)

## --------------------------------------------------------------------------
print(annotated_result)

## --------------------------------------------------------------------------
print(head(as.data.frame(annotated_result)))

## --------------------------------------------------------------------------
# Use preloaded tfbs from package sample_data. Could be manually loaded as with:
# tfbs_file = system.file('extdata','tfbs.bed.gz', package = 'methylSig')
# tfbs = rtracklayer::import(tfbs_file, genome = 'hg19')

print(tfbs)

## --------------------------------------------------------------------------
# Significance threshold
dmcList = result$fdr < 0.05 & abs(result$meth.diff) > 25

# Perform the test
tfbs_enrichment = methylSig.tfbsEnrichTest(myDiff = result, dmcList = dmcList, tfbsInfo = tfbs)

# Take a look at the first few rows
print(head(tfbs_enrichment))

