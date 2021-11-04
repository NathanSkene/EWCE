# First, check the number of annotations equals the number of columns
# in the expression data
check_annotLevels <- function(annotLevels,
    exp) {
    err_msg2 <- paste0(
        "Error: length of all annotation levels must equal",
        " the number of columns in exp matrix"
    )
    out <- lapply(annotLevels, test <- function(x, exp) {
        if (length(x) != dim(exp)[2]) {
            stop(err_msg2)
        }
    }, exp)
}
