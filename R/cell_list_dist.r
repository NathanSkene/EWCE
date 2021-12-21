#' cell_list_dist
#' 
#' specificity is generated in the main_CellTypeAnalysis_Preperation.r file
#' 
#' @param hitGenes List of gene symbols containing the target gene list. 
#' @inheritParams bootstrap_enrichment_test
#' @return Output (#fix).
#' 
#' @keywords internal
cell_list_dist <- function(hitGenes,
                           sct_data,
                           annotLevel) {
    ValidGenes <-
        rownames(sct_data[[annotLevel]]$specificity)[
            rownames(sct_data[[annotLevel]]$specificity) %in% hitGenes
        ]
    temp <- sct_data[[annotLevel]]$specificity[ValidGenes, ]

    # If the function was based a single gene list... just return temp
    # if(is.null(dim(hitGenes)[1])){
    #    return(temp)
    # }else{
    return(apply(temp, 2, sum))
    # }
}
