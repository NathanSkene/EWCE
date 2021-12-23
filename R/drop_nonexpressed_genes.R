#' Drop genes with zero counts
#'
#' Remove rows (genes) in which counts sum to zero.
#'
#' @param exp Gene expression matrix.
#' @param verbose Print messages.
#' 
#' @return List of filtered \code{exp}.
#'
#' @keywords internal
#' @importFrom Matrix rowSums
drop_nonexpressed_genes <- function(exp,
                                    verbose = TRUE) {
    messager("Checking for non-expressed genes.", v = verbose)
    orig.dims <- dim(exp)
    row.sums <- Matrix::rowSums(exp) # MUST be from Matrix
    n_zeros <- sum(row.sums <= 0, na.rm = TRUE)
    #### Drop genes ####
    if (n_zeros > 0) {
        exp <- exp[row.sums > 0, ]
        messager(nrow(exp) - orig.dims[1], "/", nrow(exp),
            "non-expressed genes dropped",
            v = verbose
        )
    }
    return(exp)
}
