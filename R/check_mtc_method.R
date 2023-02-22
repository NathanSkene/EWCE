check_mtc_method <- function(mtc_method){
    err_msg <- paste0(
        "ERROR: Invalid mtc_method argument. Please see",
        " '?p.adjust' for valid methods."
    )
    if (!mtc_method %in% c(
        stats::p.adjust.methods
    )) {
        stop(err_msg)
    }
}