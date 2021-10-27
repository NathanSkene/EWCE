create_quadrants <- function(data_byGene2) {
    #### GET QUANTILES FOR TRANSCRIPT LENGTH + GC CONTENT ####
    tl_quants <- stats::quantile(data_byGene2$transcript_length,
        probs = seq(0.1, 1, 0.1)
    )
    gc_quants <- stats::quantile(data_byGene2$percentage_gene_gc_content,
        probs = seq(0.1, 1, 0.1)
    )
    #### ASSIGN EACH GENE TO A QUANTILE QUADRANT ####
    quadrant <- matrix(0,
        nrow = dim(data_byGene2)[1],
        ncol = 2
    )
    colnames(quadrant) <- c("TL", "GC")
    for (i in seq_len(dim(data_byGene2)[1])) {
        quadrant[i, 1] <- which(data_byGene2[i, 2] < tl_quants)[1]
        quadrant[i, 2] <- which(data_byGene2[i, 3] < gc_quants)[1]
    }
    data_byGene2$uniq_quad <- sprintf("%s_%s", quadrant[, 1], quadrant[, 2])
    data_byGene2 <- data_byGene2[
        !data_byGene2$uniq_quad %in% c("2_NA", "NA_2", "3_NA"),
    ]
    return(data_byGene2)
}
