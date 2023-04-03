#' cell_list_dist
#' 
#' specificity is generated in the main_CellTypeAnalysis_Preperation.r file
#' 
#' @param hits List of gene symbols containing the target gene list. 
#' @inheritParams bootstrap_enrichment_test
#' @returns The summed specificity of each celltype
#'  across a  set of \code{hits}.
#' 
#' @keywords internal
cell_list_dist <- function(hits,
                           sct_data,
                           annotLevel) {
    ValidGenes <-
        rownames(sct_data[[annotLevel]]$specificity)[
            rownames(sct_data[[annotLevel]]$specificity) %in% hits
        ]
    temp <- sct_data[[annotLevel]]$specificity[ValidGenes, ,drop=FALSE]

    # If the function was based a single gene list... just return temp
    # if(is.null(dim(hits)[1])){
    #    return(temp)
    # }else{
    return(apply(temp, 2, sum))
    # }
}
