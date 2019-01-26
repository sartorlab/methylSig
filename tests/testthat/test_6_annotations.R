utils::data(sample_data, package = 'methylSig')

test_that('Test annotations', {
    dmcList = msig_cpgs$fdr < 0.05 & abs(msig_cpgs$meth.diff) > 25
    myDiff_annotated = methylSigAnnotation(myDiff = msig_cpgs, dmcList = dmcList, annotations = cpg_annots)

    expect_match(class(myDiff_annotated), 'GRanges')
    expect_true(all(genome(myDiff_annotated) == 'hg19'))
})
