calculate_specificity_for_level <- function(ctd_oneLevel,
                                            matrix_name = "mean_exp",
                                            as_sparse = TRUE,
                                            verbose = TRUE) {
    if (!matrix_name %in% names(ctd_oneLevel)) {
        messager(matrix_name, "not found in ctd_oneLevel.", v = verbose)
    }
    expMatrix <- ctd_oneLevel[[matrix_name]]
    normalised_meanExp <- Matrix::t(Matrix::t(expMatrix) *
        (1 / colSums(expMatrix)))
    spec <- normalised_meanExp /
        (apply(normalised_meanExp, 1, sum, na.rm = TRUE) + 1e-12)
    ctd_oneLevel$specificity <- to_sparse_matrix(
        exp = spec,
        as_sparse = as_sparse,
        verbose = verbose
    )
    return(ctd_oneLevel)
}
