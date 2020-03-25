#' Group cytosine / CpG level data into regions based on genomic regions
#'
#' An optional function to aggregate cytosine / CpG level data into regions based on a \code{GRanges} set of genomic regions.
#'
#' @param bs a \code{BSseq} object.
#' @param group_column a \code{character} string indicating the column of \code{pData(bs)} to use for determining group membership.
#' @param min_samples_per_group a named \code{integer} vector indicating the minimum number of samples with non-zero coverage required for maintaining a locus.
#'
#' @return A \code{BSseq} object with only those loci having \code{min_samples_per_group}.
#'
#' @examples
#' data(BS.cancer.ex, package = 'bsseqData')
#'
#' filter_loci_by_group_coverage(
#'     bs = BS.cancer.ex,
#'     group_column = 'Type',
#'     min_samples_per_group = c('cancer' = 3, 'normal' = 3)
#' )
#'
#' @export
filter_loci_by_group_coverage = function(bs, group_column, min_samples_per_group) {

    # Check missing
    if (missing(bs)) {
        stop('Must pass bs as a BSseq object.')
    }
    if (missing(group_column)) {
        stop('Must pass group_column as a character string.')
    }
    if (missing(min_samples_per_group)) {
        stop('Must pass min_samples_per_group as a named integer vector.')
    }

    #####################################

    # Check types
    if (!is(bs, 'BSseq')) {
        stop('bs must be class BSseq.')
    }
    if (!(is(group_column, 'character') && length(group_column) == 1)) {
        stop('group_column must be a character string.')
    }
    if (!is(min_samples_per_group, 'numeric')) {
        stop('min_samples_per_group must be a named integer vector.')
    }

    #####################################

    # Check valid group_column name
    if (!(group_column %in% colnames(pData(bs)))) {
        stop(sprintf('group_column: %s not in column names of pData(bs): %s',
            group_column, paste(colnames(pData(bs)), collapse = ', ')))
    }

    # Check valid factor names in group_column of pData(bs)
    if (!(all(names(min_samples_per_group) %in% pData(bs)[, group_column]))) {
        stop(sprintf('Not all names of min_samples_per_group are in group_column: %s',
            paste(setdiff(names(min_samples_per_group), pData(bs)[, group_column]), collapse = ', ') ))
    }

    #####################################

    # Extract sample names belonging to each group given. NOTE, min_sample_per_group
    # allows users to give more than two groups, and will require all group
    # minimums are satisfied. Likely, it will most often be the case that
    # there will only be two groups
    group_samples = lapply(names(min_samples_per_group), function(f){
        pData(bs)[, group_column] == f
    })
    names(group_samples) = names(min_samples_per_group)

    # Previous filter_ functions have set sample/loci coverages not meeting
    # the filtering requirements to 0, use this fact.
    # NOTE, coercion to DelayedArray::DelayedArray because it seems that for
    # tiled data (or smaller data), a DelayedArray isn't returned by getCoverage()
    logical_cov_mat = DelayedArray::DelayedArray(bsseq::getCoverage(bs, type = 'Cov') > 0)

    # rowSums of the logical matrix subsetted on the group columns should equal
    # or exceed the corresponding group's min_samples_per_group
    keep_group_loci = lapply(names(group_samples), function(group) {
        DelayedMatrixStats::rowSums2(
            x = logical_cov_mat,
            cols = group_samples[[group]],
            value = TRUE, na.rm = TRUE) >= min_samples_per_group[group]
    })
    names(keep_group_loci) = names(min_samples_per_group)

    # Entry-wise and() and will give loci matching all group thresholds
    keep_loci = Reduce(`&`, keep_group_loci)

    # Check that there are some loci to keep. Say which groups were too strict.
    if(!any(keep_loci)) {
        zero_groups = vapply(keep_group_loci, sum, 1, USE.NAMES = TRUE) == 0

        stop(sprintf('Thresholds for the following groups were too strict: %s.
            Relax thresholds for these groups and try filtering again.',
            paste(names(zero_groups), collapse = ', ')))
    }

    return(bs[keep_loci])
}
