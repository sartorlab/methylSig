#' Plot DM status distribution in annotations
#'
#' Generates a bar plot showing the distribution of differential methylation status of loci in selected annotations.
#'
#' @param myAnnots A \code{GRanges} object resulting from running \code{methylSigAnnotation} on the result of \code{methylSigCalc}.
#' @param annotation_order A character vector which orders and subsets the annotations for the plot.
#' @param status_order A character vector which orders and subsets the DM status.
#' @param position One of \code{fill}, \code{stack}, \code{dodge}. See \code{ggplot2} documentation for details.
#' @param plot_title The plot label.
#' @param legend_title The legend label.
#' @param x_label The x-axis label.
#' @param y_label The y-axis label.
#'
#' @return A \code{ggplot2} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' # Annotate the msig_cpgs results
#' data(data, package = 'methylSig')
#'
#' # Use the genome of msig_cpgs and build annotations for CpG features
#' genome = GenomeInfoDb::genome(msig_cpgs)
#' annots = annotatr::build_annotations(genome = genome, annotations = paste(genome, c('cpgs'), sep='_'))
#'
#' # Decide what counts as differentially methylated
#' dmcList = msig_cpgs$fdr < 0.05 & abs(msig_cpgs$meth.diff) > 25
#'
#' # Annotate
#' myDiff_annotated = methylSigAnnotation(myDiff = msig_cpgs, dmcList = dmcList, annotations = annots)
#'
#' # Set the order vectors
#' cpg_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter')
#' dm_order = c('DR','DS','No DM')
#'
#' methylSigPlotStatus(myAnnots = myDiff_annotated, annotation_order = cpg_order, status_order = dm_order,
#'     position = 'fill', plot_title = 'DM Status in CpG Annots.', legend_title = 'Annotations',
#'     x_label = 'DM Status', y_label = 'Proportion')
#'
#' @export
methylSigPlotStatus = function(myAnnots, annotation_order = NULL, status_order = NULL, position = 'fill',
    plot_title = 'DM Status', legend_title = 'Annotations', x_label = 'DM Status', y_label = 'Proportion') {

    plot = annotatr::plot_categorical(
        annotated_regions = myAnnots,
        x = 'dm_status',
        fill = 'annot.type',
        x_order = status_order,
        fill_order = annotation_order,
        position = position,
        plot_title = plot_title,
        legend_title = legend_title,
        x_label = x_label,
        y_label = y_label)

    return(plot)
}

#' Data visualization function
#'
#' Generates data visualization plot of methylation data for a specified genomic interval.
#'
#'   This function offers a unique two-tiered visualization of the methylation data depending on the zoom level. For narrow regions (<1mbp) where at most 500 CpG sites have data reads, users can visualize sample-specific coverage levels and percent methylation at each site, together with group averages, significance levels and a number of genomic annotations.
#'
#' @param myAnnots A \code{GRanges} object resulting from running \code{methylSigAnnotation} on the result of \code{methylSigCalc}.
#' @param annotation_order A character vector which orders and subsets the annotations for the plot.
#' @param status_order A character vector which orders and subsets the DM status.
#' @param bin_width A vector of two numbers (from, to) to specify the region to visualize on chromosome `chr'.
#' @param plot_title The plot label.
#' @param legend_title The legend label.
#' @param x_label The x-axis label.
#' @param y_label The y-axis label.
#' @param legend_facet_label Label explaining the meaning of the gray bars in the resulting histogram.
#' @param legend_cum_label Label explaining the meaning of the red outlined bars in the resulting histogram.
#'
#' @return A \code{ggplot2} object which can be viewed by calling it, or saved with \code{ggplot2::ggsave}.
#'
#' @examples
#' # Annotate the msig_cpgs results
#' data(data, package = 'methylSig')
#'
#' # Use the genome of msig_cpgs and build annotations for CpG features
#' genome = GenomeInfoDb::genome(msig_cpgs)
#' annots = annotatr::build_annotations(genome = genome, annotations = paste(genome, c('cpgs'), sep='_'))
#'
#' # Decide what counts as differentially methylated
#' dmcList = msig_cpgs$fdr < 0.05 & abs(msig_cpgs$meth.diff) > 25
#'
#' # Annotate
#' myDiff_annotated = methylSigAnnotation(myDiff = msig_cpgs, dmcList = dmcList, annotations = annots)
#'
#' # Set the order vectors
#' cpg_order = c('hg19_cpg_islands','hg19_cpg_shores','hg19_cpg_shelves','hg19_cpg_inter')
#' status_order = c('hyper','hypo','none')
#'
#' methylSigPlotDiff(myAnnots = myDiff_annotated, annotation_order = cpg_order, status_order = status_order,
#'     bin_width = 10, plot_title = 'Meth. Diff. in CpG Annots.', x_label = 'DM Status', y_label = 'Proportion',
#'     legend_facet_label = 'Meth. Diff. in Annotation', legend_cum_label = 'Meth. Diff. Overall')
#'
#' @export
methylSigPlotDiff = function(myAnnots, annotation_order = NULL, status_order = NULL, bin_width = 10,
    plot_title = 'Methylation Differences', x_label = 'Methylation Difference', y_label = 'Density',
    legend_facet_label = 'Meth. Diff. in Annotation', legend_cum_label = 'Meth. Diff. Overall') {

    facet_order = list(annotation_order, status_order)

    plot = annotatr::plot_numerical(
        annotated_regions = myAnnots,
        x = 'meth.diff',
        facet = c('annot.type', 'dm_status'),
        facet_order = facet_order,
        bin_width = bin_width,
        plot_title = plot_title,
        x_label = x_label,
        y_label = y_label,
        legend_facet_label = legend_facet_label,
        legend_cum_label = legend_cum_label)

    return(plot)
}
