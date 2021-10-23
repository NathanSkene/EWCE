#' Run DGE: \pkg{limma}
#'
#' Run Differential Gene Expression with \pkg{limma}.
#'
#' @return \code{limma} results
#' @inheritParams drop_uninformative_genes
#'
#' @keywords internal
#' @importFrom limma lmFit eBayes
#' @importFrom stats model.matrix
run_limma <- function(exp,
                      level2annot,
                      verbose = TRUE,
                      ...) {
    messager("DGE:: Limma...", v = verbose)
    ## Prepare groupings
    level2_options <- as.factor(as.character(level2annot))
    mod_matrix <- stats::model.matrix(~level2_options)
    fit <- limma::lmFit(exp, mod_matrix, ...)
    eb <- limma::eBayes(fit)
    return(eb)
}
