messager <- function(..., v = TRUE) {
    if (v) {
        msg <- paste(...)
        message(msg)
    }
}
