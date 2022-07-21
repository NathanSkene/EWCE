#' Get CTD matrix names
#' 
#' Find the names of all data matrices in a CellTypeDataset.
#' 
#' @param ctd CellTypeDataset. If set to \code{NULL} (default), 
#' will simply return all possible matrix names. 
#' @param matrices Matrix names to search for.
#' @param verbose Print messages.
#' @keywords internal
#' @returns List of matrix names.
get_ctd_matrix_names <- function(ctd=NULL,
                                 matrices=c("mean_exp", 
                                            "median_exp",
                                            "specificity", 
                                            "median_specificity", 
                                            "specificity_quantiles"),
                                 verbose=TRUE){ 
    if(is.null(ctd)){
        messager("Returning all possible matrix names.",v=verbose)
        return(matrices)
    }
    nms <- lapply(ctd, function(ctd_lvl){
        names(ctd_lvl)[names(ctd_lvl) %in% matrices]
    }) |> unlist() |> unique()
    messager("Found",length(nms),"matrix types across",
             length(ctd),"CTD levels.",v=verbose)
    return(nms)
}
