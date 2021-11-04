#' bin_specificity_into_quantiles
#'
#' \code{bin_specificity_into_quantiles} is an internal function used to convert
#' add '$specificity_quantiles' to a ctd
#'
#' @param ctdIN A single annotLevel of a ctd, i.e. ctd[[1]] (the function is
#' intended to be used via apply)
#' @param numberOfBins Number of quantile 'bins' to use (40 is recommended)
#' @param as_sparse Convert to sparseMatrix.
#' @param verbose Print messages.
#'
#' @return A ctd with $specificity_quantiles.
#' @examples
#' library(ewceData)
#' ctd <- ctd()
#' ctd <- lapply(ctd, bin_specificity_into_quantiles, numberOfBins = 40)
#' print(ctd[[1]]$specificity_quantiles[1:3, ])
#' @export
#' @importFrom methods as
bin_specificity_into_quantiles <- function(ctdIN,
    numberOfBins,
    as_sparse = TRUE,
    verbose = TRUE) {
    specQ <- apply(ctdIN$specificity, 2,
        FUN = bin_columns_into_quantiles,
        numberOfBins = numberOfBins
    )
    ctdIN$specificity_quantiles <- to_sparse_matrix(
        exp = specQ,
        as_sparse = as_sparse,
        verbose = verbose
    )
    rownames(ctdIN$specificity_quantiles) <- rownames(ctdIN$specificity)
    return(ctdIN)
}
