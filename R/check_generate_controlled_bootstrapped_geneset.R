#' generate_controlled_bootstrap_geneset
#' 
#' Check input arguments to \link[EWCE]{generate_controlled_bootstrap_geneset}.
#' 
#' @inheritParams generate_controlled_bootstrap_geneset
#' @inheritParams bootstrap_enrichment_test 
#' @return Null output.
#' 
#' @keywords internal
check_generate_controlled_bootstrap_geneset <- function(controlledCT,
                                                        annotLevel,
                                                        sct_data,
                                                        hitGenes) {
    err_msg <- paste0(
        "ERROR: controlledCT cannot be NULL in",
        " generate_controlled_bootstrap_geneset"
    )
    if (is.null(controlledCT)) {
        stop(err_msg)
    }
    err_msg2 <- paste0(
        "ERROR: annotLevel cannot be greater than the number",
        " of annotation levels in sct_data"
    )
    if (annotLevel > length(sct_data)) {
        stop(err_msg2)
    }
    # Check all controlledCT are in single cell data
    err_msg3 <- paste0(
        "ERROR: not all controlledCT are in",
        " colnames(sct_data[[annotLevel]]$specificity)"
    )
    if (sum(!controlledCT %in%
        colnames(sct_data[[annotLevel]]$specificity)) != 0) {
        stop(err_msg3)
    }
    err_msg4 <- paste0(
        "ERROR: length(hitGenes)==0. Perhaps your gene list is",
        " from the wrong species? It should be converted to",
        " orthologs of the same species as the single cell",
        " dataset."
    )
    if (length(hitGenes) == 0) {
        stop(err_msg4)
    }
}
