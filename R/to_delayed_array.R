#' Convert object to DelayedArray
#'
#' Convert a variety of object types to DelayedArray format.
#'
#' @param exp Object.
#' @param as_DelayedArray Whether to convert \code{exp} to DelayedArray.
#' @param verbose Print messages.
#'
#' @keywords internal
#' @importFrom DelayedArray DelayedArray
to_delayed_array <- function(exp,
                             as_DelayedArray = TRUE,
                             verbose = TRUE) {
    if (as_DelayedArray && (!is_delayed_array(exp))) {
        messager("Converting to DelayedArray.", v = verbose)
        if (!is_matrix(exp)) {
            exp <- as.matrix(exp)
        }
        exp <- DelayedArray::DelayedArray(exp)
    }
    return(exp)
}
