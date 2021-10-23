#' Prepare dendrogram
#'
#' \code{prep_dendro} adds a dendrogram to a CellTypeDataset (CTD).
#'
#' @param ctdIN A single annotLevel of a ctd, i.e. ctd[[1]] (the function is
#' intended to be used via apply).
#' @return A CellTypeDataset with dendrogram plotting info added.
#'
#' @keywords internal
#' @importFrom stats hclust dist
#' @importFrom Matrix t
prep_dendro <- function(ctdIN) {
    requireNamespace("ggdendro")
    # euclidean distances between the rows
    binned_file_dist <- stats::dist(Matrix::t(ctdIN$specificity_quantiles))
    binned_file_dist_hclust <- stats::hclust(binned_file_dist)
    ddata <- ggdendro::dendro_data(binned_file_dist_hclust,
        type = "rectangle"
    )
    ordered_cells <- as.character(ddata$labels$label)
    a1 <- ggplot(ggdendro::segment(ddata)) +
        geom_segment(aes_string(
            x = "x", y = "y",
            xend = "xend", yend = "yend"
        )) +
        coord_flip() +
        ggdendro::theme_dendro()
    a1 <- a1 + scale_x_continuous(expand = c(0, 1.3))
    b1 <- ggplot(ggdendro::segment(ddata)) +
        geom_segment(aes_string(
            x = "x", y = "y",
            xend = "xend", yend = "yend"
        )) +
        ggdendro::theme_dendro()
    b1 <- b1 + scale_x_continuous(expand = c(0, 1.3))
    ctdIN$plotting <- list()
    ctdIN$plotting$ggdendro_vertical <- a1
    ctdIN$plotting$ggdendro_horizontal <- b1
    ctdIN$plotting$cell_ordering <- ordered_cells
    return(ctdIN)
}

#' prep.dendro
#'
#' @inherit prep_dendro
prep.dendro <- function(ctdIN) {
    .Deprecated("prep_dendro")
    ctdIN <- prep_dendro(ctdIN = ctdIN)
    return(ctdIN)
}
