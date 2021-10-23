filter_variance_quantiles <- function(exp,
                                      level2annot,
                                      n_quantiles = 10,
                                      min_variance_quantile = 5,
                                      verbose = TRUE) {
    messager("+ Filtering by variance quantiles...", v = verbose)
    #### Calculate gene variance across cell types means ####
    gene_variance <- stats::setNames(
        DelayedMatrixStats::rowVars(exp),
        rownames(exp)
    )
    #### Convert to qunatiles ####
    quant <- calc_quantiles(
        v = gene_variance,
        n_quantiles = n_quantiles,
        report_filters = verbose
    )
    #### Remove genes below the min_variance_quantile ####
    gene_variance <- gene_variance[quant >= (min_variance_quantile / n_quantiles)]
    # DelayedMatrixStats::rowQuantiles(gene_variance)
    # hist(log(gene_variance), breaks=50)
    orig.dim <- dim(exp)
    messager(paste(
        nrow(exp) - length(gene_variance), "/", nrow(exp),
        "genes dropped @ DGE min_variance_qunatile >=", min_variance_quantile
    ), v = verbose)
    exp <- exp[names(gene_variance), ]
    return(exp)
}
