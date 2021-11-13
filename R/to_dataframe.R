#' Convert object to data.frame
#'
#' Convert a variety of object types to data.frame format.
#'
#' @param X Object.
#' @param verbose Print messages.
#'
#' @keywords internal
#' @importFrom methods is
to_dataframe <- function(X,
                         verbose = TRUE) {
    if (methods::is(X, "data.frame")) {
        return(X)
    } else if (is_matrix(X)) {
        messager("Converting to data.frame", v = verbose)
        nn <- colnames(X)
        rr <- rownames(X)
        if (is_delayed_array(X) || is_sparse_matrix(X)) {
            X <- as.matrix(X)
        }
        X <- data.frame(X,
            stringsAsFactors = FALSE,
            check.rows = FALSE,
            check.names = FALSE
        )
        colnames(X) <- nn
        rownames(X) <- rr
        return(X)
    } else {
        stop("Format ", methods::is(X)[1], " is not supported.")
    }
}
