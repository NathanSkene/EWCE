#' Calculate quantiles
#' 
#' @param vec Numeric vector.
#' @param n_quantiles Number of quantile bins to use.
#' Defaults to deciles (\code{n_quantiles=10}).
#' @param verbose Print messages.
#'
#' @return Quantiles.
#'
#' @keywords internal
#' @importFrom stats setNames ecdf
calc_quantiles <- function(vec,
                            n_quantiles = 10,
                            verbose = TRUE) {
    calc_percentile <- stats::ecdf(seq(1, n_quantiles))
    quantiles <- calc_percentile(vec)
    # Report the number of items left after filtering at each quantile
    if (verbose) {
        messager("N remaining at each quantile filter:")
        genes_left <- lapply(unique(quantiles), function(x) {
            stats::setNames(length(vec[quantiles >= x]), paste0(">=", x*10))
        }) %>%
            unlist() %>%
            sort() %>%
            rev()
        print(genes_left)
    }
    return(quantiles)
}
