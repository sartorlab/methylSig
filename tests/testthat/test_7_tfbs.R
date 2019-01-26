utils::data(sample_data, package = 'methylSig')

test_that('Test annotations', {
    dmcList = msig_cpgs$fdr < 0.05 & abs(msig_cpgs$meth.diff) > 25
    tfbs_results = methylSig.tfbsEnrichTest(myDiff = msig_cpgs, dmcList = dmcList, tfbsInfo = tfbs)

    expect_equal(sum(tfbs_results$pvalue < 0.05), 4)
})
