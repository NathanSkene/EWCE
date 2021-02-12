#' convert_new_ewce_to_old
#'
#' \code{convert_new_ewce_to_old} Used to get an old style EWCE ctd file 
#' from a new one
#'
#' @param ctd A celltype data structure containing $mean_exp and $specificity
#' @param lvl The annotation level to extract
#' @return celltype_data The old style data structure
convert_new_ewce_to_old <- function(ctd, lvl) {
    celltype_data <- list()
    celltype_data[[1]] <- list()
    celltype_data[[1]]$cell_dists <- ctd[[lvl]]$specificity
    celltype_data[[1]]$all_scts <- ctd[[lvl]]$mean_exp
    return(celltype_data)
}
