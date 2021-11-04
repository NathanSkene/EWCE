#' calculate_meanexp_for_level
#'
#' @return One level of a CellTypeDataset.
#'
#' @keywords internal
#' @importFrom stats model.matrix
calculate_meanexp_for_level <- function(ctd_oneLevel,
                                        expMatrix,
                                        as_sparse = TRUE,
                                        verbose = TRUE) {
    err_msg <- paste0(
        "There are an equal number of cell types in expMatrix",
        " and ctd_oneLevel but the names do not match"
    )
    if (dim(expMatrix)[2] == length(unique(ctd_oneLevel$annot))) {
        message(dim(expMatrix)[2])
        message(length(ctd_oneLevel$annot))
        if (sum(!colnames(expMatrix) == ctd_oneLevel$annot) != 0) {
            stop(err_msg)
        }
        ctd_oneLevel$mean_exp <- to_sparse_matrix(
            exp = expMatrix,
            as_sparse = as_sparse,
            verbose = verbose
        )
    } else {
        # Sum reads in each cell type
        mm <- stats::model.matrix(~ 0 + ctd_oneLevel$annot)
        colnames(mm) <- names(table(ctd_oneLevel$annot))
        mat.summary.mm1 <- expMatrix %*% mm

        # Divide by the number of cells to get the mean
        cellCounts <- table(ctd_oneLevel$annot)
        for (i in seq_len(dim(mat.summary.mm1)[2])) {
            mat.summary.mm1[, i] <- mat.summary.mm1[, i] / cellCounts[i]
        }
        ctd_oneLevel$mean_exp <- to_sparse_matrix(
            exp = mat.summary.mm1,
            as_sparse = as_sparse,
            verbose = verbose
        )
    }
    return(ctd_oneLevel)
}
