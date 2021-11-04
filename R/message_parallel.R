#' Print messages
#'
#' Print messages even from within parallelised functions.
#'
#' @param ... Message input.
#'
#' @keywords internal
message_parallel <- function(...) {
    system(sprintf('echo "%s"', paste0(..., collapse = "")))
}
