#' \code{bin_columns_into_quantiles}
#'
#' \code{bin_columns_into_quantiles} is an internal function used to convert a
#' vector of specificity  into a vector of specificity quantiles.
#' This function can be iterated across a matrix using \link[base]{apply} 
#' to create a matrix of specificity quantiles.
#' @param vec The vector of gene of specificity values.
#' @param numberOfBins Number of quantile bins to use (40 is recommended).
#' @param defaultBin Which bin to assign when there's only one
#' non-zero quantile. In situations where there's only one non-zero quantile,
#' \link[base]{cut} throws an error. Avoid these situations by
#' using a default quantile.
#' @returns A vector with same length as \code{vec} but with columns storing
#' quantiles instead of specificity.
#' @examples
#' ctd <- ewceData::ctd()
#' ctd[[1]]$specificity_quantiles <- apply(ctd[[1]]$specificity, 2,
#'     FUN = bin_columns_into_quantiles)
#' @export
#' @importFrom stats quantile
bin_columns_into_quantiles <- function(vec,
                                       numberOfBins = 40,
                                       defaultBin = as.integer(
                                           numberOfBins / 2)
                                       ) {
    
    quantileValues <- rep(0, length(vec))
    breaks <- unique(stats::quantile(vec[vec > 0],
        probs = seq(0, 1, by = 1 / numberOfBins),
        na.rm = TRUE
    )) 
    if (length(breaks) > 1) {
        quantileValues[vec > 0] <- as.numeric(cut(vec[vec > 0],
            breaks = breaks,
            include.lowest = TRUE
        ))
    } else {
        ## In situations where there's only one non-zero quantile, 
        ##  cut() throws an error.
        ## Avoid these situations by using a default quantile.
        messager(
            "+ <2 non-zero quantile bins detected in column.",
            "Assigning these values to default quantile ",
            "(", defaultBin, ")"
        )
        quantileValues[vec > 0] <- defaultBin
    }
    return(quantileValues)
}
