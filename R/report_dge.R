#' Report DGE
#' 
#' Report differential gene expression (DGE) results
#' 
#' @param exp Gene expression matrix.
#' @param keep_genes Genes kept after DGE.
#' @inheritParams drop_uninformative_genes
#' 
#' @keywords internal
report_dge <- function(exp, 
                       keep_genes, 
                       adj_pval_thresh = .05,
                       verbose = TRUE){
    messager(paste(
        formatC(nrow(exp) - length(keep_genes), big.mark = ","),
        "/",
        formatC(nrow(exp), big.mark = ","),
        "genes dropped @ DGE adj_pval_thresh <", adj_pval_thresh
    ), v = verbose)
}
