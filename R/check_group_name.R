#' Check group name
#' 
#' Ensure \code{groupName} argument is provided to 
#' \link[EWCE]{generate_celltype_data}.
#' @inheritParams generate_celltype_data
#' @return Null output.
#' 
#' @keywords internal
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
