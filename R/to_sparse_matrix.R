#' Convert object to sparse matrix
#'
#' Convert a variety of object types to sparse matrix format.
#'
#' @param exp Object.
#' @param as_sparse Whether to convert \code{exp} to sparse matrix
#' @param verbose Print messages.
#' 
#' @keywords internal
#' @importFrom DelayedArray DelayedArray
to_sparse_matrix <- function(exp,
                             as_sparse = TRUE,
                             verbose = TRUE) {
    if (as_sparse) {
        messager("Converting to sparse matrix.", v = verbose)
        if (!is_sparse_matrix(exp)) {
            if (!is_matrix(exp)) {
                exp <- as.matrix(exp)
            }
            exp <- methods::as(exp, "sparseMatrix")
        }
    } else {
        #### Convert to dense matrix ####
        exp <- as.matrix(exp, "matrix")
        #### Convert characters to numbers ####
        exp <- check_numeric(exp = exp)
    }
    return(exp)
}
