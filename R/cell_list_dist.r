# specificity is generated in the main_CellTypeAnalysis_Preperation.r file
cell_list_dist <- function(hitGenes, sct_data, annotLevel) {
    ValidGenes <- 
        rownames(sct_data[[annotLevel]]$specificity)[
                rownames(sct_data[[annotLevel]]$specificity) %in% hitGenes]
    temp <- sct_data[[annotLevel]]$specificity[ValidGenes, ]

    # If the function was based a single gene list... just return temp
    # if(is.null(dim(hitGenes)[1])){
    #    return(temp)
    # }else{
    return(apply(temp, 2, sum))
    # }
}
