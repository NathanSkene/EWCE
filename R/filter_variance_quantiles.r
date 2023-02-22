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
#' 
#' @returns Filtered \code{exp}.
#' @keywords internal
filter_variance_quantiles <- function(exp,
                                      log10_norm = TRUE,
                                      n_quantiles = 10,
                                      min_variance_quantile = as.integer(
                                          n_quantiles / 2
                                      ),
                                      verbose = TRUE) {
    # templateR:::args2vars(filter_variance_quantiles)
    # exp <- ewceData::ctd()[[1]]$mean_exp
    
    exp_orig <- exp
    messager("Filtering by variance quantiles.", v = verbose)
    #### Log normalise to avoid skewed quantiles ####
    if (isTRUE(log10_norm)) {
        exp <- log10(exp + 1e-12)
    }
    #### Convert to DelayedArray to take advantage of rowVars func #####
    exp <- to_delayed_array(exp, verbose = verbose)
    #### Calculate gene variance across cell types means ####
    gene_variance <- stats::setNames(
        DelayedMatrixStats::rowVars(exp),
        rownames(exp)
    )
    #### Convert to quantiles ####
    quant <- bin_columns_into_quantiles(vec = gene_variance,
                                        numberOfBins = n_quantiles)
    #### Remove genes below the min_variance_quantile ####
    gene_variance <- gene_variance[quant>=min_variance_quantile] 
    messager(paste(
        formatC(nrow(exp) - length(gene_variance), big.mark = ","),
        "/",
        formatC(nrow(exp), big.mark = ","),
        "genes dropped @ DGE min_variance_quantile >=", min_variance_quantile
    ), v = verbose)
    #### Return filtered original data ####
    return(exp_orig[names(gene_variance), ])
}
