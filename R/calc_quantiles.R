#' Calculate quantiles
#'
#' @return Quantiles.
#'
#' @keywords internal
#' @importFrom stats setNames
calc_quantiles <- function(v,
                           n_quantiles = 10,
                           report_filters = TRUE) {
    calc_percentile <- stats::ecdf(seq(1, n_quantiles))
    quantiles <- calc_percentile(v)
    # Report the number of items left after filtering at each quantile
    if (report_filters) {
        messager("N remaining at each quantile filter:")
        genes_left <- lapply(unique(quantiles), function(x) {
            stats::setNames(length(v[quantiles >= x]), paste0(">=", x))
        }) %>%
            unlist() %>%
            sort() %>%
            rev()
        print(genes_left)
    }
    return(quantiles)
}
