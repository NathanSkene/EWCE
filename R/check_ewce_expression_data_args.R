check_ewce_expression_data_args <- function(sortBy,
                                            tt,
                                            thresh){
    err_msg <- paste0(
        "ERROR: tt does not contain a column with value",
        " passed in sortBy argument"
    )
    # Check the arguments
    if (!sortBy %in% colnames(tt)) {
        stop(err_msg)
    }
    err_msg2 <- paste0(
        "ERROR: length of table is less than",
        " twice the size of threshold"
    )
    if (dim(tt)[1] < (thresh * 2)) {
        stop(err_msg2)
    } 
}