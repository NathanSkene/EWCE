#' Filter genes in a CellTypeDataset
#'
#' Removes rows from each matrix within a CellTypeDataset (CTD) that are not 
#' within \code{gene_subset}.
#'
#' @param ctd CellTypeDataset.
#' @param gene_subset Genes to subset to.
#'
#' @returns Filtered CellTypeDataset.
#'
#' @export  
#' @examples  
#' ctd <- ewceData::ctd()
#' ctd <- standardise_ctd(ctd, input_species="mouse")
#' gene_subset <- rownames(ctd[[1]]$mean_exp)[1:100]
#' ctd_subset <- EWCE::filter_ctd_genes(ctd = ctd, gene_subset = gene_subset) 
filter_ctd_genes <- function(ctd,
                             gene_subset) {
    message("Filtering CTD to ", 
            formatC(length(gene_subset),big.mark = ","), " genes.")
    
    new_ctd <- lapply(get_ctd_levels(ctd), 
                      function(lvl) {
        message("level: ", lvl)
        mat_names <- get_ctd_matrix_names(ctd[lvl])                  
        ctd_lvl <- ctd[[lvl]] 
        other_names <- names(ctd_lvl)[!names(ctd_lvl) %in% mat_names]
        new_ctd_lvl <- lapply(mat_names, function(mat_nm) {
            message("   - ", mat_nm)
            ctd_lvl[[mat_nm]][rownames(ctd_lvl[[mat_nm]]) %in% gene_subset, ]
        }) |> `names<-`(mat_names)
        for (nm in other_names) {
            new_ctd_lvl[nm] <- ctd_lvl[nm]
        }
        return(new_ctd_lvl)
    })
    return(new_ctd)
}
