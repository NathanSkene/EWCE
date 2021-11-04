#' Assess whether an object is a DelayedArray.
#'
#' Assess whether an object is a DelayedArray or one of
#'  its derived object types.
#'
#' @param X Object.
#' @importFrom methods is
is_delayed_array <- function(X) {
    methods::is(X, "DelayedMatrix") |
        methods::is(X, "DelayedArray")
}
