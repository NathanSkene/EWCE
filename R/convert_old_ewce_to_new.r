#' convert_old_ewce_to_new
#'
#' \code{convert_old_ewce_to_new} Used to get an new style EWCE ctd file
#' (mean_exp/specificity) from old ones (all_scts).
#'
#' If you've already loaded it and want to pass it as a celltype_data
#' structure, then don't set level1 or level2.
#'
#' @param level1 File path to old level1 of EWCE ctd.
#' @param level2 File path to old level2 of EWCE ctd.
#' @param celltype_data The ctd to be converted.
#' 
#' @return CellTypeData in the new data structure style.
#' 
#' @keywords internal
convert_old_ewce_to_new <- function(level1 = NA,
                                    level2 = NA,
                                    celltype_data = NA) {
    ctd <- list()
    if (!is.na(level1)) {
        celltype_data <- load_rdata(level1)

        ctd[[1]] <- list()
        ctd[[1]]$specificity <- celltype_data[[1]]$cell_dists
        ctd[[1]]$mean_exp <- celltype_data[[1]]$all_scts
        if("annot" %in% names(celltype_data[[1]])){
            ctd[[1]]$annot <- celltype_data[[1]]$annot
        }
        
        if(!is.na(level2)){
            celltype_data <- load_rdata(level2)
            ctd[[2]] <- list()
            ctd[[2]]$specificity <- celltype_data[[1]]$cell_dists
            ctd[[2]]$mean_exp <- celltype_data[[1]]$all_scts
            if("annot" %in% names(celltype_data[[1]])){
                ctd[[2]]$annot <- celltype_data[[1]]$annot
            }
        } 
    } else {
        for (i in seq_len(length(celltype_data))) {
            ctd[[i]] <- list()
            ctd[[i]]$specificity <- celltype_data[[i]]$cell_dists
            ctd[[i]]$mean_exp <- celltype_data[[i]]$all_scts
            if("annot" %in% names(celltype_data[[i]])){
                ctd[[i]]$annot <- celltype_data[[i]]$annot
            }
        }
    }
    return(ctd)
}
