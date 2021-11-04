#' Filter variance quantiles
#'
#' Remove rows in \code{exp} that do not vary substantially across rows.
#'
#' @param exp Gene expression matrix.
#' @param log10_norm Log10-normalise \code{exp} before computing variance.
#' @param n_quantiles Number of quantile bins to use.
#' Defaults to deciles (\code{n_quantiles=10}).
#' @param min_variance_quantile The minimum variance quantile
#'  to keep values from.
#' @param verbose Print messages.
filter_variance_quantiles <- function(exp,
                                      log10_norm = TRUE,
                                      n_quantiles = 10, 
                                      min_variance_quantile = as.integer(
                                          n_quantiles / 2
                                      ),
                                      verbose = TRUE) {
    exp_orig <- exp
    messager("Filtering by variance quantiles.", v = verbose)
    #### Log normalise to avoid skewed quantiles ####
    if(log10_norm){
        exp <- log10(exp + 1e-12)
    }
    #### Convert to DelayedArray to take advantage of rowVars func #####
    exp <- to_delayed_array(exp, verbose = verbose)
    #### Calculate gene variance across cell types means ####
    gene_variance <- stats::setNames(DelayedMatrixStats::rowVars(exp),
                                     rownames(exp))
    #### Convert to quantiles ####
    quant <- calc_quantiles(
        vec = gene_variance,
        n_quantiles = n_quantiles,
        verbose = verbose
    )
    #### Remove genes below the min_variance_quantile ####
    gene_variance <- gene_variance[
        quant >= (min_variance_quantile / n_quantiles)
    ]
    # DelayedMatrixStats::rowQuantiles(gene_variance)
    # hist(log(gene_variance), breaks=50)
    orig.dim <- dim(exp)
    messager(paste(
        formatC(nrow(exp) - length(gene_variance), big.mark = ","),
        "/",
        formatC(nrow(exp),big.mark = ","),
        "genes dropped @ DGE min_variance_quantile >=", min_variance_quantile
    ), v = verbose)
    #### Return filtered original data ####
    return(exp_orig[names(gene_variance), ])
}
