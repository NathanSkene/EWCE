#' Assess whether an object is a sparse matrix
#'
#' Assess whether an object is a sparse matrix or one of
#'  its derived object types.
#'
#' @param X Object.
#' @importFrom methods is
is_sparse_matrix <- function(X) {
    methods::is(X, "sparseMatrix") |
        methods::is(X, "dgCMatrix") |
        methods::is(X, "dgRMatrix") |
        methods::is(X, "dgTMatrix") |
        methods::is(X, "dgeMatrix") |
        methods::is(X, "lgCMatrix")
}
