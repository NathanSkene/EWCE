check_group_name <- function(groupName) {
    err_msg3 <- paste0(
        "ERROR: groupName must be set. groupName is used to",
        " label the files created by this function."
    )
    if (is.null(groupName)) {
        stop(err_msg3)
    }
    if (is.null(groupName) || groupName == "") {
        stop(err_msg3)
    }
}
