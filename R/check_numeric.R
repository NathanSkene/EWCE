#' Check numeric
#' 
#' Ensure that a matrix is numeric. If not, it will be converted to numeric.
#' @param exp Input matrix.
#' @return Numeric expression matrix.
#' 
#' @keywords internal
#' @importFrom methods is
check_numeric <- function(exp) {
    if (methods::is(exp[1, 1], "character")) {
        storage.mode(exp) <- "numeric"
    }
    return(exp)
}
