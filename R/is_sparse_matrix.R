is_sparse_matrix <- function(X) {
    methods::is(X, "sparseMatrix") |
        methods::is(X, "dgCMatrix") |
        methods::is(X, "dgRMatrix") |
        methods::is(X, "dgTMatrix") |
        methods::is(X, "dgeMatrix") |
        methods::is(X, "lgCMatrix")
}
