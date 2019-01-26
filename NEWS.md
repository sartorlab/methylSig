# CHANGES IN VERSION 0.5.2

    - Reintroduce Seqinfo via the `assembly` parameter of `methylSigReadData`
        - If `NA`, allow user to continue, but warn that tiling and annotation cannot work unless the user manually assigns `Seqinfo` to the `BSseq` object.
        - If `assembly` is unsupported in `GenomeInfoDb::fetchExtendedChromInfoFromUCSC`, act as if it was set to `NA`.
        - If `assembly` is supported, use it.
        - `seqinfo` persists through tiling, tests for different methylation, and annotation.
    - Update built-in data based on `Seqinfo` fixes.

# CHANGES IN VERSION 0.5.1

    - Minor bugfix in methylSigTile() function. Cannot reproduce everywhere, but sometimes when a chromosome has no entries in the meth object, the tiling function failed. Now check for entries on the chromosome.

# CHANGES IN VERSION 0.5.0

    NOTE: This version of methylSig only works for R <= 3.4.4 and Bioc <= 3.6. A feature change in the bsseq Bioconductor package in Bioc >= 3.7 does not allow BSseq-class objects whose GRanges are not points, and this breaks the tiling functionality of methylSig.

    ## NEW FEATURES

    - Annotations are now done with the annotatr Bioconductor package.

    ## USER-LEVEL CHANGES

    - methylSig v0.5.0 reuses Bioconductor classes rather than the home-spun classes of earlier versions. This will improve maintainability greatly.
    - The methylSigReadData() function now is a wrapper for the bsseq::read.bismark() function, obviating the need to transform the input data in anyway. The output is a BSseq-class object.
    - As before, filtering for common SNPs (hg19 only), minCount, and maxCount are available. Destranding also remains.
    - The result of any of the tests for differential methylation are now GenomicRanges-class objects.
    - Built-in example data is now known as sample_data.

    ## BUG FIXES

    - Fixed a mistake in methylSig.tfbsEnrichTest() that mistakenly referred referred to tfbsInfo parameter as tfbs.

    ## REMOVED FEATURES

    - Removed plotting functions for retooling.
    - Do not use seqinfo for any objects. Instead, when tiling a genome, find the maximum length of each existing chromosome, add 1000, and use that as the input to GenomicRanges::tileGenome().
