#' bin.specificity.into.quantiles
#'
#' \code{bin.specificity.into.quantiles} is an internal function used to convert
#' add '$specificity_quantiles' to a ctd
#'
#' @param ctdIN A single annotLevel of a ctd, i.e. ctd[[1]] (the function is 
#' intended to be used via apply)
#' @param numberOfBins Number of quantile 'bins' to use (40 is recommended)
#' @return A ctd with $specificity_quantiles
#' @examples
#' ctd <- ctd()
#' ctd <- lapply(ctd, bin.specificity.into.quantiles, numberOfBins = 40)
#' print(ctd[[1]]$specificity_quantiles[1:3, ])
#' @export
bin.specificity.into.quantiles <- function(ctdIN, numberOfBins) {
    ctdIN$specificity_quantiles <- apply(ctdIN$specificity, 2, 
                                            FUN = bin.columns.into.quantiles, 
                                            numberOfBins = numberOfBins)
    rownames(ctdIN$specificity_quantiles) <- rownames(ctdIN$specificity)
    return(ctdIN)
}
