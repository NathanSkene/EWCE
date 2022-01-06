#' Prepare differential gene expression table
#' 
#' Prepare differential gene expression table for 
#' \link[EWCE]{generate_bootstrap_plots_for_transcriptome} or
#' \link[EWCE]{ewce_expression_data}.
#' 
#' @inheritParams ewce_expression_data
#' @inheritParams orthogene::convert_orthologs
#' 
#' @returns List of 3 items
#' 
#' @keywords internal
prepare_tt <- function(tt,
                       tt_genecol = NULL,
                       ttSpecies,
                       output_species,
                       method = "homologene",
                       verbose = TRUE){
    #### Infer gene column in tt ####
    if(is.null(tt_genecol)){
        tt_genecol <- colnames(tt)[1]   
        messager("Using 1st column of tt as gene column:",tt_genecol,
                 v=verbose)
    }
    #### Convert orthologs if needede ####
    if(ttSpecies!=output_species){
        tt <- orthogene::convert_orthologs(gene_df = tt,
                                           gene_input = tt_genecol,
                                           gene_output = "columns",
                                           input_species = ttSpecies,
                                           output_species = output_species,
                                           method = method, 
                                           verbose = verbose)
        tt_genecol <- "ortholog_gene"
        ttSpecies <- output_species
    }
    return(list(tt = tt,
                tt_genecol = tt_genecol,
                ttSpecies = ttSpecies))
}