#' Check whether a CellTypeDataset is standardised
#'
#' Check whether a CellTypeDataset was previously standardised
#' using \link[EWCE]{standardise_ctd}.
#'
#' @param ctd CellTypeDataset.
#'
#' @return Whether the \code{ctd} is standardised.
#'
#' @keywords internal
is_ctd_standardised <- function(ctd) {
    std_list <- lapply(ctd, function(x) {
        if ("standardised" %in% names(x)) {
            return(x[["standardised"]])
        } else {
            return(FALSE)
        }
    })
    is_standardised <- all(unlist(std_list) == TRUE)
    return(is_standardised)
}
