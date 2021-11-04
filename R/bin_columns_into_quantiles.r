#' \code{bin_columns_into_quantiles}
#'
#' \code{bin_columns_into_quantiles} is an internal function used to convert a
#' matrix of specificity (with columns of cell types) intom a matrix of
#' specificity quantiles
#'
#' @param matrixIn The matrix of specificity values
#' @param numberOfBins Number of quantile 'bins' to use (40 is recommended)
#' @param defaultBin Which bin to assign when there's only one
#' non-zero quantile. In situations where there's only one non-zero quantile,
#' \code{cut()} throws an error. Avoid these situations by
#' using a default quantile.
#' @return A matrix with same shape as matrixIn but with columns storing
#' quantiles instead of specificity
#' @examples
#' ctd <- ewceData::ctd()
#' ctd[[1]]$specificity_quantiles <- apply(ctd[[1]]$specificity, 2,
#'     FUN = bin_columns_into_quantiles,
#'     numberOfBins = 40
#' )
#' @export
#' @importFrom stats quantile
bin_columns_into_quantiles <- function(matrixIn,
    numberOfBins = 40,
    defaultBin = as.integer(numberOfBins / 2)) {
    quantileValues <- rep(0, length(matrixIn))
    breaks <- unique(stats::quantile(matrixIn[matrixIn > 0],
        probs = seq(0, 1, by = 1 / numberOfBins),
        na.rm = TRUE
    ))
    if (length(breaks) > 1) {
        quantileValues[matrixIn > 0] <- as.numeric(cut(matrixIn[matrixIn > 0],
            breaks = breaks,
            include.lowest = TRUE
        ))
    } else {
        ## In situations where there's only one non-zero quantile, cut() throws an error.
        ## Avoid these situations by using a default quantile.
        message(
            "+ <2 non-zero quantile bins detected in column. Assigning these values to default quantile ",
            "(", defaultBin, ")"
        )
        quantileValues[matrixIn > 0] <- defaultBin
    }
    return(quantileValues)
}
