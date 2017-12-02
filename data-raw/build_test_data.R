# In /bfx/home/rcavalca/epicore/pipeline_testing/test_errbs_1.2.3/00-methCall
# do the following in bash

# for file in *cytosine_report.txt
# do
#     file_base=`basename ${file} '_cytosine_report.txt'`
#     echo ${file_base}
#     grep chr21 ${file} | head -200000 | gzip > ${file_base}_chr21_cytosine_report.txt.gz
# done

devtools::load_all()

# Get data
files = list.files(
    path = '~/ActiveProjects/EpiCore/pipeline_testing/test_errbs_1.2.3/00-methCall',
    pattern = 'txt.gz', full.names = TRUE)

sample.ids = basename(files)
sample.ids = gsub('_chr21_cytosine_report.txt.gz', '', sample.ids)

pData = data.frame(
    Sample_Names = sample.ids,
    DR_vs_DS = relevel(factor(c('DR','DS','DR','DS','DR','DS')), ref = 'DS'),
    row.names = sample.ids,
    stringsAsFactors = FALSE)

data = methylSigReadData(
    fileList = files,
    pData = pData,
    assembly = 'hg19',
    destranded = TRUE,
    maxCount = 500,
    minCount = 10,
    filterSNPs = TRUE,
    num.cores = 4,
    fileType = 'cytosineReport')

save(data, file = 'data/data.RData', compress = 'xz')
