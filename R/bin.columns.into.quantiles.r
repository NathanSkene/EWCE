#' bin.columns.into.quantiles
#'
#' \code{bin.columns.into.quantiles} is an internal function used to convert a 
#' matrix of specificity (with columns of cell types) intom a matrix of 
#' specificity quantiles
#'
#' @param matrixIn The matrix of specificity values
#' @param numberOfBins Number of quantile 'bins' to use (40 is recommended)
#' @return A matrix with same shape as matrixIn but with columns storing 
#' quantiles instead of specificity
#' @examples
#' data(ctd, package="ewceData")
#' ctd[[1]]$specificity_quantiles <- apply(ctd[[1]]$specificity, 2,
#'     FUN = bin.columns.into.quantiles,
#'     numberOfBins = 40
#' )
#' @export
bin.columns.into.quantiles <- function(matrixIn, numberOfBins = 40) {
    quantileValues <- rep(0, length(matrixIn))
    quantileValues[matrixIn > 0] <- as.numeric(cut(matrixIn[matrixIn > 0],
        breaks = unique(quantile(matrixIn[matrixIn > 0], 
                                    probs = seq(0, 1, by = 1 / numberOfBins), 
                                    na.rm = TRUE)),
        include.lowest = TRUE
    ))
    return(quantileValues)
}
