#' prep.dendro
#'
#' \code{prep.dendro} adds a dendrogram to a cts
#'
#' @param ctdIN A single annotLevel of a ctd, i.e. ctd[[1]] (the function is 
#' intended to be used via apply)
#' @return A ctd with dendrogram plotting info added
#' @examples
#' ctd <- ctd()
#' ctd <- lapply(ctd, EWCE::bin.specificity.into.quantiles, numberOfBins = 40)
#' ctd <- lapply(ctd, EWCE::prep.dendro)
#' @export
#' @import ggdendro
prep.dendro <- function(ctdIN) {
    # euclidean distances between the rows
    binned_file_dist <- dist(t(ctdIN$specificity_quantiles)) 
    binned_file_dist_hclust <- hclust(binned_file_dist)
    ddata <- ggdendro::dendro_data(binned_file_dist_hclust, type = "rectangle")
    ordered_cells <- as.character(ddata$labels$label)
    a1 <- ggplot(segment(ddata)) +
        geom_segment(aes_string(x = "x", y = "y", 
                                    xend = "xend", yend = "yend")) +
        coord_flip() +
        theme_dendro()
    a1 <- a1 + scale_x_continuous(expand = c(0, 1.3))
    b1 <- ggplot(segment(ddata)) +
        geom_segment(aes_string(x = "x", y = "y", 
                                    xend = "xend", yend = "yend")) +
        theme_dendro()
    b1 <- b1 + scale_x_continuous(expand = c(0, 1.3))
    ctdIN$plotting <- list()
    ctdIN$plotting$ggdendro_vertical <- a1
    ctdIN$plotting$ggdendro_horizontal <- b1
    ctdIN$plotting$cell_ordering <- ordered_cells
    return(ctdIN)
}
