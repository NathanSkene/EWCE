#' Drop cells with zero gene counts
#'
#' Remove columns (cells) in which (gene) counts sum to zero.
#'
#' @param exp Gene expression matrix.
#' @param annotLevels Cell-wise annotations to be subset
#' if some cells are dropped.
#' @param verbose Print messages.
#'
#' @keywords internal
#' @importFrom Matrix colSums
drop_nonexpressed_cells <- function(exp,
    annotLevels,
    verbose = TRUE) {
    messager("Checking for cells with no expressed genes.", v = verbose)
    orig.dims <- dim(exp)
    col.sums <- Matrix::colSums(exp) # MUST be from Matrix
    n_zeros <- sum(col.sums <= 0, na.rm = TRUE)
    #### Drop genes ####
    if (n_zeros > 0) {
        exp <- exp[, col.sums > 0]
        ### Make sure annotLevels are also subset appropriately.
        annotLevels <- lapply(annotLevels, function(x) x[col.sums > 0])
        messager(nrow(exp) - orig.dims[1], "/", nrow(exp),
            "cells dropped",
            v = verbose
        )
    }
    return(list(
        exp = exp,
        annotLevels = annotLevels
    ))
}
