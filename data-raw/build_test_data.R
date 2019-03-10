# In /bfx/home/rcavalca/epicore/pipeline_testing/test_errbs_1.2.3/00-methCall
# do the following in bash

# for file in *cytosine_report.txt
# do
#     file_base=`basename ${file} '_cytosine_report.txt'`
#     echo ${file_base}
#     grep chr21 ${file} | head -200000 | gzip > ${file_base}_chr21_cytosine_report.txt.gz
# done
#
# for file in `ls *cytosine_report*`
# do
#     file_base=`basename ${file} '_chr21_cytosine_report.txt.gz'`
#     echo ${file_base}
#     gunzip -c ${file} | awk -v OFS='\t' '$4 + $5 > 5 {print $0}' | gzip > ${file_base}.txt.gz
# done

devtools::load_all()

# Get data
files = list.files(
    path = 'inst/extdata',
    pattern = 'MDAMB', full.names = TRUE)

sample.ids = basename(files)
sample.ids = gsub('.txt.gz', '', sample.ids)

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
    num.cores = 2)

tiled_meth = methylSigTile(meth = meth, tiles = NULL, win.size = 1000)

msig_cpgs = methylSigCalc(
    meth = meth,
    comparison = 'DR_vs_DS',
    dispersion = 'both',
    local.info = FALSE,
    local.winsize = 200,
    min.per.group = c(3,3),
    weightFunc = methylSig_weightFunc,
    T.approx = TRUE,
    num.cores = 1)

msig_tiles = methylSigCalc(
    meth = tiled_meth,
    comparison = 'DR_vs_DS',
    dispersion = 'both',
    local.info = FALSE,
    local.winsize = 200,
    min.per.group = c(3,3),
    weightFunc = methylSig_weightFunc,
    T.approx = TRUE,
    num.cores = 1)

tfbs_file = system.file('extdata','tfbs.bed.gz', package = 'methylSig')
tfbs = rtracklayer::import(tfbs_file, genome = 'hg19')

cpg_annots = annotatr::build_annotations(genome = 'hg19', annotations = 'hg19_cpg_islands')

save(meth, tiled_meth, msig_cpgs, msig_tiles, tfbs, cpg_annots, file = 'data/sample_data.RData', compress = 'xz')
