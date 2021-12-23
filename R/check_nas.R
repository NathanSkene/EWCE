#' Check NAs
#' 
#' Check for any NAs in an expression matrix.
#' 
#' @param exp Expression matrix.
#' @return Null output.
#' 
#' @keywords internal
check_nas <- function(exp) {
    err_msg <- paste0(
        "NA values detected in expresson matrix. All NA values",
        " should be removed before running EWCE."
    )
    if (sum(is.na(exp)) > 0) {
        stop(err_msg)
    }
}
