calculate_specificity_for_level <- function(ctd_oneLevel,
                                            as_sparse = TRUE,
                                            verbose = TRUE) {
    #### Guards against issues with DelayedArray
    expMatrix <- as.matrix(ctd_oneLevel$mean_exp)
    normalised_meanExp <- Matrix::t(Matrix::t(expMatrix) *
        (1 / colSums(expMatrix)))
    spec <- normalised_meanExp /
        (apply(normalised_meanExp, 1, sum) + 0.000000000001)
    ctd_oneLevel$specificity <- to_sparse_matrix(
        exp = spec,
        as_sparse = as_sparse,
        verbose = verbose
    )
    return(ctd_oneLevel)
}
