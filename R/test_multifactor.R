#' Calculates differential methylation statistics under general experimental design
#'
#' @param meth A \code{BSseq-class} object to calculate differential methylation statistics. See \code{methylSigReadData} for how to read in methylation data.
#' @param design A \code{data.frame} for experimental design. Should contain as many rows as there are columns (samples) in \code{meth}.
#' @param formula A formula for the linear model. It should refer to column names from \code{design}. NOTE: The intercept is included by default if omitted. One can omit the intercept with a formula such as \code{'~ 0 + group'}. For clarity, it helps to include the intercept explicitly as in \code{'~ 1 + group'}.
#' @param contrast A contrast matrix for hypothesis testing. The number of rows should match the number of columns \code{design}.
#' @param group.term A string indicating which term in \code{formula} contains group information on which to apply the \code{min.per.group} parameter. Currently assumes that this factor contains ONLY TWO LEVELS.
#' @param min.per.group A vector with two numbers specifying the minimum number of samples required to perform the test for differential methylation. If it is a single number, both groups will use it as the minimum requried number of samples. Default is \code{c(3,3)}.  NOTE: The ordering of this parameter with respect to the groups should be \code{c(reference, other)}, where \code{reference} refers to the reference level in the \code{design[, group.term]} factor.
#'
#' @return A \code{GRanges} object containing the following \code{mcols}:
#' \describe{
#'   \item{stat}{ The dispersion estimate. }
#'   \item{pvalue}{ The log likelihood ratio. }
#'   \item{fdr}{ Degrees of freedom used when \code{T.approx = TRUE}. }
#'   \item{mean_methylation_columns}{ Mean methylation for each factor in each column of design, when there are fewer than 5 factors in the factor. }
#' }
#'
#' @seealso \code{\link{methylSigReadData}}
#'
#' @examples
#' data(data, package = 'methylSig')
#'
#' # Example with implicit intercept
#' design1 = data.frame(group = bsseq::pData(data)$DR_vs_DS)
#' contrast1 = matrix(c(0,1), ncol = 1)
#' result1 = methylSigDSS(
#'     meth = data,
#'     design = design1,
#'     formula = '~ group',
#'     contrast = contrast1,
#'     group.term = 'group',
#'     min.per.group=c(3,3))
#'
#' # Example with no intercept and corresponding change in contrast
#' # NOTE: result2 is the same as result1
#' contrast2 = matrix(c(-1,1), ncol = 1)
#' result2 = methylSigDSS(
#'     meth = data,
#'     design = design1,
#'     formula = '~ 0 + group',
#'     contrast = contrast2,
#'     group.term = 'group',
#'     min.per.group=c(3,3))
#'
#' # Example with subject pairing
#' design2 = data.frame(
#'     group = bsseq::pData(data)$DR_vs_DS,
#'     subject = factor(c(1,1,2,2,3,3)))
#' contrast2 = matrix(c(0,1,0,0), ncol = 1)
#' result2 = methylSigDSS(
#'     meth = data,
#'     design = design2,
#'     formula = '~ group + subject',
#'     contrast = contrast2,
#'     group.term = 'group',
#'     min.per.group=c(3,3))
#'
#' @keywords differentialMethylation
#'
#' @export
methylSigDSS = function(
    meth,
    design,
    formula,
    contrast,
    group.term,
    min.per.group=c(3,3)) {

    #####################################

    # Determine the formula components so we only return mean methylation information
    # for items in the formula intersected with criteria below
    # TO DO: Add support for interaction terms
    formula_components = unlist(strsplit(formula, '[+]'))
    formula_components = gsub(' ', '', formula_components)
    formula_components = gsub('~', '', formula_components)
    formula_components = formula_components[formula_components != '0']

    ############################################################################
    # Filter according to min.per.group parameter
    if(length(min.per.group) == 1) {
        min.per.group = c(min.per.group,min.per.group)
    }

    #####################################
    # Get the group labels, NOTE: THIS ASSUMES CORRECT REFERENCE LEVEL SET
    group2 = levels(design[, group.term])[2]
    group1 = levels(design[, group.term])[1]

    # Determine which rows of pData belong to which group
    # / which columns of Cov and M matrices belong to which group
    group2_idx = which(design[, group.term] == group2)
    group1_idx = which(design[, group.term] == group1)

    #####################################
    # Determine which sites are valid to test according to min.per.group
    all_cov = as.matrix(bsseq::getCoverage(meth, type = 'Cov'))

    # Determine which loci satisfy min.per.group
    valid_idx = which(
        base::rowSums(all_cov[, group2_idx] > 0) >= min.per.group[2] & base::rowSums(all_cov[, group1_idx] > 0) >= min.per.group[1]
    )

    meth = meth[valid_idx]

    ############################################################################
    # Test for differential methylation

    dss_fit = DSS::DMLfit.multiFactor(BSobj = meth, design = design, formula = as.formula(formula), smoothing=FALSE)

    # Need to remove any rows with NAs because the contrast will fail
    # design; fit: beta, var.beta; formula; gr; X
    na_idx = which(is.na(dss_fit$fit$beta[,1]))
    if(length(na_idx) > 0) {
        dss_fit$fit$beta = dss_fit$fit$beta[-na_idx, ]
        dss_fit$fit$var.beta = dss_fit$fit$var.beta[-na_idx, ]
        dss_fit$gr = dss_fit$gr[-na_idx]
    }

    test_result = DSS::DMLtest.multiFactor(DMLfit = dss_fit, Contrast = contrast)

    # Create the GRanges return object and harmonize column names with methylSigCalc()
    dss_result = dss_fit$gr
    GenomicRanges::mcols(dss_result) = test_result[,c('stat','pvals','fdrs')]
    colnames(GenomicRanges::mcols(dss_result)) = c('stat','pvalue','fdr')

    # Retain metadata from dss_fit, design, formula, and contrast
    dss_result_metadata = list(
        design = design,
        formula = formula,
        contrast = contrast,
        beta_fit = dss_fit$fit$beta,
        var_beta_fit = dss_fit$fit$var.beta,
        X = dss_fit$X
    )
    S4Vectors::metadata(dss_result) = dss_result_metadata

    ############################################################################
    # Recover group mean methylation

    # Construct methylation means matrix to add as mcols to dss_result
    all_meth = as.matrix(bsseq::getCoverage(meth, type = 'M'))
    all_cov = as.matrix(bsseq::getCoverage(meth, type = 'Cov'))
    perc_meth = all_meth / all_cov
    if(length(na_idx) > 0) {
        perc_meth = perc_meth[-na_idx, ]
    }

    #####################################
    # Collect the correct columns of dss_fit$design to use for column groups

    # Want the column to be a factor
    factor_cols = sapply(dss_fit$design, class) == 'factor'
    # Want the column to have 5 or fewer levels (if not factor returns FALSE)
    level_cols = sapply(lapply(dss_fit$design, levels), length) <= 5
    # Which column names are valid (above two criteria plus in formula)?
    char_idx = intersect(names(which(factor_cols & level_cols)), formula_components)

    #####################################

    # Determine the column indices for each factor level in char_idx
    col_idxs = lapply(char_idx, function(idx){
        col = dss_fit$design[, idx]
        groups = unique(as.character(col))
        tmp = lapply(groups, function(group){
            which(col == group)
        })
        names(tmp) = groups

        return(tmp)
    })
    names(col_idxs) = char_idx

    #####################################

    # Take the row means over the selected columns
    row_means = lapply(col_idxs, function(group){
        sapply(group, function(idxs){
            base::rowMeans(perc_meth[, idxs], na.rm = TRUE)
        })
    })

    means_df = Reduce(cbind, row_means)
    colnames(means_df) = paste('mean_meth', colnames(means_df), sep= '.')

    #####################################

    # Add meth difference based on group.term

    means_df = cbind(
        meth.diff = means_df[, paste('mean_meth', group2, sep='.')] - means_df[, paste('mean_meth', group1, sep='.')],
        means_df)

    #####################################

    # Add the mean methylations to the dss_result
    GenomicRanges::mcols(dss_result) = cbind(GenomicRanges::mcols(dss_result), means_df)

    ############################################################################

    return(dss_result)
}
