#' Extract significant results
#'
#' Extract significant results from output of
#' \link[EWCE]{bootstrap_enrichment_test}.
#'
#' @param full_results Output of \link[EWCE]{bootstrap_enrichment_test}.
#' @param q_threshold Maximum multiple-testing-corrected p-value to include.
#' @inheritParams bootstrap_enrichment_test
#'
#' @keywords internal
#' @importFrom stats p.adjust
get_sig_results <- function(full_results,
                            mtc_method = "BH",
                            q_threshold = .05,
                            verbose = TRUE) {
    res <- full_results$results
    if (!"q" %in% colnames(res)) {
        res$q <- stats::p.adjust(res$p, method = mtc_method)
    }
    messager(nrow(res), "signficiant enrichment results @", mtc_method, "<",
        q_threshold,
        v = verbose
    )
    return(res)
}
