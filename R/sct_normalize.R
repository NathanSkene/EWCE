#' Normalize expression matrix
#'
#' Normalize expression matrix by accounting for library size.
#' Uses \pkg{sctransform}.
#'
#' @param exp Gene x cell expression matrix.
#' @param as_sparse Convert \code{exp} to sparse matrix.
#' @param verbose Print messages.
#'
#' @return Normalised expression matrix.
#'
#' @examples
#' cortex_mrna <- ewceData::cortex_mrna()
#' exp_sct_normed <- EWCE::sct_normalize(exp = cortex_mrna$exp[1:300, ])
#' @export
#' @importFrom Matrix t colSums
sct_normalize <- function(exp,
                          as_sparse = TRUE,
                          verbose = TRUE) {
    requireNamespace("sctransform")
    exp <- to_sparse_matrix(
        exp = exp,
        as_sparse = as_sparse,
        verbose = verbose
    )
    sct <- sctransform::vst(
        umi = exp,
        return_cell_attr = TRUE,
        verbosity = if (verbose) 2 else 0
    )
    exp_sct <- sctransform::correct_counts(
        x = sct,
        umi = exp,
        # UMI_corrected
        verbosity = if (verbose) 2 else 0
    )
    exp_sct_normed <- Matrix::t(Matrix::t(exp_sct) *
        (1 / Matrix::colSums(exp_sct)))
    return(exp_sct_normed)
}
