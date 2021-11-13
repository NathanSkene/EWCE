check_full_results <- function(full_results,
                               sct_data) {
    if (!is.null(full_results)) {
        err_msg <- paste0(
            "ERROR: full_results is not valid output from the",
            " bootstrap_enrichment_test function"
        )
        if (length(full_results) != 3) {
            stop(err_msg)
        }
        err_msg2 <- paste0(
            "ERROR: No cell types in full_results are found in",
            " sct_data. Perhaps the wrong annotLevel was used?"
        )
        if (sum(!as.character(unique(full_results$results$CellType)) %in%
            colnames(sct_data[[1]]$specificity)) ==
            length(as.character(unique(full_results$results$CellType)))) {
            stop(err_msg2)
        }
        if (sum(!as.character(unique(full_results$results$CellType)) %in%
            colnames(sct_data[[1]]$specificity)) > 0) {
            stop(err_msg2)
        }
    }
}
