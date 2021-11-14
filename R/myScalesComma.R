#' \code{myScalesComma}
#'
#' Adjusts \pkg{ggplot2} label display. See \link[scales]{comma} for details.
#' Support function for \link[EWCE]{plot_log_bootstrap_distributions}.
#' 
#' @return Numeric vector
#' 
#' @keywords internal
#' @importFrom scales comma
myScalesComma <- function(x) {
    requireNamespace("scales")
    return(scales::comma(x = x, accuracy = 0.01))
}
