#' Check whether object is a CellTypeDataset
#'
#' Check whether an object is a CellTypeDataset.
#'
#' @param ctd Object.
#'
#' @return boolean
#' 
#' @keywords internal
is_celltypedataset <- function(ctd) {
    (!is.function(ctd)) &&
        all(c("annot", "mean_exp", "specificity") %in% names(ctd[[1]]))
}
