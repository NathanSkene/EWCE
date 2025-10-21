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
prep_dendro <- function(ctdIN, 
                        expand=c(0, 0.66)) {
    requireNamespace("ggplot2")
    requireNamespace("ggdendro")
    # euclidean distances between the rows
    binned_file_dist <- stats::dist(Matrix::t(ctdIN$specificity_quantiles))
    binned_file_dist_hclust <- stats::hclust(binned_file_dist)
    ddata <- ggdendro::dendro_data(binned_file_dist_hclust,
        type = "rectangle"
    )
    ordered_cells <- as.character(ddata$labels$label)
    #### Vertical dendrogram ####
    a1 <- ggplot2::ggplot(ggdendro::segment(ddata)) +
        ggplot2::geom_segment(
            ggplot2::aes(
            x = .data$x, y = .data$y,
            xend = .data$xend, yend = .data$yend
        )) +
        ggplot2::coord_flip() +
        ggdendro::theme_dendro()
    if(!is.null(expand)){
        a1 <- a1 + ggplot2::scale_x_continuous(expand = expand)
    }
    #### Horizontal dendrogram ####
    b1 <- ggplot(ggdendro::segment(ddata)) +
        ggplot2::geom_segment(ggplot2::aes(
            x = .data$x, y = .data$y,
            xend = .data$xend, yend = .data$yend
        )) +
        ggdendro::theme_dendro()
    if(!is.null(expand)){
        b1 <- b1 + ggplot2::scale_x_continuous(expand = expand)    
    }
    #### Make nested list ####
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
