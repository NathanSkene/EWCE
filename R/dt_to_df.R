dt_to_df <- function(exp) {
    if (methods::is(exp, "data.table")) {
        col1 <- colnames(exp)[1]
        exp <- data.frame(exp,
            row.names = col1,
            check.rows = FALSE,
            check.names = FALSE
        )
    }
    return(exp)
}
