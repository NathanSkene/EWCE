check_numeric <- function(exp) {
    if (methods::is(exp[1, 1], "character")) {
        storage.mode(exp) <- "numeric"
    }
    return(exp)
}
