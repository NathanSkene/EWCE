is_celltypedataset <- function(ctd) {
    (!is.function(ctd)) &&
        all(c("annot", "mean_exp", "specificity") %in% names(ctd[[1]]))
}
