get_ctd_levels <- function(ctd,
                           max_only = FALSE) {
    # This is necessary in case further meta-data such as $name is used
    if (!is.null(names(ctd))) {
        lvls <- names(ctd)
    } else {
        lvls <- seq(1, length(ctd))
    }
    if (max_only) {
        return(max(lvls))
    } else {
        return(lvls)
    }
}
