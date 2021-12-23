#' Assess whether an object is a Matrix
#'
#' Assess whether an object is a Matrix or one of
#'  its derived object types.
#'
#' @param X Object.
#' 
#' @return boolean
#' 
#' @importFrom methods is
is_matrix <- function(X) {
    methods::is(X, "Matrix") || methods::is(X, "matrix")
}
