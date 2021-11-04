#' convert_new_ewce_to_old
#'
#' \code{convert_new_ewce_to_old} Used to get an old style EWCE ctd file
#' from a new one
#'
#' @param ctd A cell type data structure containing
#'  "mean_exp" and "specificity".
#' @param lvl The annotation level to extract.
#' 
#' @return CellTypeData in the old data structure style.
#' 
#' @keywords internal
convert_new_ewce_to_old <- function(ctd,
                                    lvl) {
    celltype_data <- list()
    celltype_data[[1]] <- list()
    celltype_data[[1]]$cell_dists <- ctd[[lvl]]$specificity
    celltype_data[[1]]$all_scts <- ctd[[lvl]]$mean_exp
    celltype_data[[1]]$annot <- ctd[[lvl]]$annot
    return(celltype_data)
}
