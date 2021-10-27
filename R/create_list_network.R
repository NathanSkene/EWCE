create_list_network <- function(data_byGene2,
                                hitGenes_NEW,
                                numBOOT = 10000) {
    # Get all sctSpecies genes in each quadrant
    quad_genes <- list()
    for (uq in unique(data_byGene2$uniq_quad)) {
        quad_genes[[uq]] <-
            unique(data_byGene2[data_byGene2$uniq_quad == uq, "HGNC.symbol"])
    }

    list_network <- matrix("",
        nrow = numBOOT,
        ncol = length(hitGenes_NEW)
    )
    count <- 0
    for (gene in hitGenes_NEW) {
        count <- count + 1
        this_gene_quad <- data_byGene2[
            data_byGene2$HGNC.symbol == gene,
            "uniq_quad"
        ][1]
        candidates <- as.vector(unlist(quad_genes[this_gene_quad]))
        list_network[, count] <- sample(
            x = candidates,
            size = numBOOT,
            replace = TRUE
        )
    }
    return(list_network)
}
