get_sig_results <- function(full_results,
                            mtc_method = "BH",
                            q_threshold = .05,
                            verbose = TRUE) {
    res <- full_results$results
    res$q <- stats::p.adjust(res$p, method = mtc_method)
    messager(nrow(res), "signficiant enrichment results @", mtc_method, "<", q_threshold, v = verbose)
    return(res)
}
