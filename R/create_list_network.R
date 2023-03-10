#' \code{create_list_network}
#'
#' Support function for \code{prepare_genesize_control_network}.
#'
#' @return List network
#'
#' @keywords internal
#' @importFrom  parallel mclapply
#' @importFrom data.table data.table
create_list_network <- function(data_byGene2,
                                hits_NEW,
                                reps = 10000,
                                no_cores = 1) {

    # Get all sctSpecies genes in each quadrant
    quad_genes <- list()
    for (uq in unique(data_byGene2$uniq_quad)) {
        quad_genes[[uq]] <-
            unique(data_byGene2[data_byGene2$uniq_quad == uq, "HGNC.symbol"])
    }

    list_network <- parallel::mclapply(hits_NEW, function(gene) {
        this_gene_quad <- data_byGene2[
            data_byGene2$HGNC.symbol == gene,
            "uniq_quad"
        ][1]
        candidates <- as.vector(unlist(quad_genes[this_gene_quad]))
        data.table::data.table(sample(
            x = candidates,
            size = reps,
            replace = TRUE
        ))
    }, mc.cores = no_cores) |>
        do.call(what = "cbind") |>
        as.matrix() |>
        `colnames<-`(NULL)
    return(list_network)
}
