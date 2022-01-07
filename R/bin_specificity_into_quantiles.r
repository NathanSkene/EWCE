#' bin_specificity_into_quantiles
#'
#' \code{bin_specificity_into_quantiles} is an internal function used to convert
#' add '$specificity_quantiles' to a ctd
#'
#' @param ctdIN A single annotLevel of a ctd, i.e. \code{ctd[[1]]} 
#' (the function is intended to be used via \code{apply}). 
#' @param numberOfBins Number of quantile 'bins' to use (40 is recommended).
#' @param matrix_name Name of the specificity matrix to create 
#' (default: "specificity_quantiles"). 
#' @param as_sparse Convert to sparseMatrix.
#' @param verbose Print messages.
#'
#' @returns A ctd with "specificity_quantiles" matrix in each level 
#' (or whatever \code{matrix_name} was set to.).
#' @examples
#' ctd <- ewceData::ctd()
#' ctd <- lapply(ctd, EWCE::bin_specificity_into_quantiles, numberOfBins = 40)
#' print(ctd[[1]]$specificity_quantiles[1:3, ])
#' @export
#' @importFrom methods as
bin_specificity_into_quantiles <- function(ctdIN,
                                           numberOfBins,
                                           matrix_name = 
                                               "specificity_quantiles",
                                           as_sparse = TRUE,
                                           verbose = TRUE) {
    specQ <- apply(ctdIN$specificity, 2,
        FUN = bin_columns_into_quantiles,
        numberOfBins = numberOfBins
    )
    ctdIN[[matrix_name]] <- to_sparse_matrix(
        exp = specQ,
        as_sparse = as_sparse,
        verbose = verbose
    )
    rownames(ctdIN[[matrix_name]]) <- rownames(ctdIN$specificity)
    return(ctdIN)
}
