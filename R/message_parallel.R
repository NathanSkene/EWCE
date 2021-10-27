message_parallel <- function(...) {
    system(sprintf('echo "%s"', paste0(..., collapse = "")))
}
